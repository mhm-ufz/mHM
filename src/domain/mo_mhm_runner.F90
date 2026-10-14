!> \file    mo_mhm_runner.F90
!> \copydoc mo_mhm_runner

!> \brief   Forward execution and startup messages for mHM v6.
!> \details Domains and their exchange state belong to each invocation. Callers
!! configure logging and message visibility before starting a forward run.
!> \version 0.1
!> \changelog
!! - Oct 2026: extracted forward execution and startup reporting from app/driver.F90.
!> \authors Sebastian Mueller
!> \date    Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_domain
#include "logging.h"
module mo_mhm_runner
  use mo_logging
  use mo_kind, only: i4
  use mo_domain, only: domain_t
  use mo_exchange_type, only: get_n_domains
  use mo_os, only: path_abspath, check_path_isdir
  use mo_utils, only: optval
  use mo_version, only: version, version_date
  implicit none
  private

  public :: run_forward, startup_message
  character(*), parameter, public :: default_main_file = "mhm.nml" !< Default main namelist shared with the CLI.
  character(*), parameter, public :: default_para_file = "parameters.nml" !< Default parameter namelist shared with the CLI.
  character(*), parameter, public :: default_out_file = "outputs.nml" !< Default output namelist shared with the CLI.

contains

  !> \brief Configure and execute a forward run from namelist paths.
  !> \details Relative paths are resolved against cwd without changing the process
  !! working directory. All domains are connected before any domain is initialized.
  subroutine run_forward(main_file, para_file, out_file, cwd)
    character(*), intent(in), optional :: main_file !< Main namelist path; defaults to mhm.nml.
    character(*), intent(in), optional :: para_file !< Parameter namelist path; defaults to parameters.nml.
    character(*), intent(in), optional :: out_file !< Output namelist path; defaults to outputs.nml.
    character(*), intent(in), optional :: cwd !< Run root; defaults to the process working directory.
    character(:), allocatable :: main_path, para_path, out_path, run_root
    character(*), parameter :: separator = repeat("-", 72)
    type(domain_t), allocatable, target :: domains(:)
    integer(i4) :: n_domains, i

    call startup_message()
    main_path = optval(main_file, default_main_file)
    para_path = optval(para_file, default_para_file)
    out_path = optval(out_file, default_out_file)
    run_root = path_abspath(optval(cwd, "."))
    call check_path_isdir(run_root, raise=.true.)

    n_domains = get_n_domains(main_file=main_path, cwd=run_root)
    log_info(*) "ALLOCATE DOMAINS: ", n_domains
    allocate(domains(n_domains))

    ! Establish every domain's configuration, restart definition, and publications first.
    do i = 1_i4, size(domains)
      log_text(*) separator
      log_info(*) "DOMAIN: ", i
      call domains(i)%create(main_file=main_path, domain=i, cwd=run_root)
      log_text(*) separator
      call domains(i)%configure(main_file=main_path, para_file=para_path, out_file=out_path)
      log_text(*) separator
      call domains(i)%prepare_restart()
      log_text(*) separator
      call domains(i)%connect()
    end do

    ! Complete and release each domain using the unchanged lifecycle and time loop.
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
      call domains(i)%destroy()
    end do
  end subroutine run_forward

  !> \brief Show the startup splash screen and OpenMP execution status.
  subroutine startup_message()
    !$ use omp_lib, only: omp_get_num_threads
    !$ integer(i4) :: n_threads
    logical :: openmp_enabled

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
    log_text(*) "|                 Lead development: Sebastian Müller                   |"
    log_text(*) "|                Originally by L. Samaniego & R. Kumar                 |"
    log_text(*) "'----------------------------------------------------------------------'"

    openmp_enabled = .false.
    !$omp parallel
    !$omp single
    !$ n_threads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    !$ openmp_enabled = .true.
    !$ log_info(*) "OpenMP enabled with ", n_threads, " threads."
    if (.not.openmp_enabled) then
      log_info(*) "OpenMP not enabled."
    end if
  end subroutine startup_message

end module mo_mhm_runner
