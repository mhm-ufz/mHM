!> \file runner.f90
!> \copydoc runner

!> \brief Forward execution for the Python package and command line.
!> \details Configures application logging and delegates the domain lifecycle to mo_mhm_runner.
!> \date Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup py_mhm
module runner
  implicit none
  private
  public :: run_forward
contains
  !> \brief Execute one forward run with the requested logging behavior.
  subroutine run_forward(main_file, para_file, out_file, cwd, verbose, quiet, &
      log_date, log_time, log_line, log_scope, log_file)
    use mo_kind, only: i4
    use mo_mhm_runner, only: execute_forward => run_forward
    use mo_logging
    use mo_message, only: use_log
    use mo_cli_support, only: set_verbosity_level
    character(*), intent(in) :: main_file !< Main namelist path relative to cwd, or absolute.
    character(*), intent(in) :: para_file !< Parameter namelist path relative to cwd, or absolute.
    character(*), intent(in) :: out_file !< Output namelist path relative to cwd, or absolute.
    character(*), intent(in) :: cwd !< Run root; does not change the process working directory.
    integer(i4), intent(in) :: verbose !< Number of increases to logging verbosity.
    integer(i4), intent(in) :: quiet !< Number of decreases to logging verbosity.
    logical, intent(in) :: log_date !< Include the current date in log messages.
    logical, intent(in) :: log_time !< Include the current time in log messages.
    logical, intent(in) :: log_line !< Include source locations in log messages.
    character(*), intent(in) :: log_scope !< Comma-separated logging scope filter.
    character(*), intent(in) :: log_file !< Log destination relative to the process cwd; empty selects the console.
    integer :: saved_unit, saved_error_unit, unit

    ! Reset relative verbosity on every invocation, including repeated Python calls.
    log_level = LOG_INFO
    call log_set_config(verbose=verbose, quiet=quiet, log_output_date=log_date, &
      log_output_time=log_time, log_output_line=log_line, log_scope_filter=log_scope)
    use_log = .true.
    call set_verbosity_level(3_i4 - quiet)
    saved_unit = log_unit
    saved_error_unit = log_unit_error
    if (len_trim(log_file) > 0) then
      open(newunit=unit, file=log_file, status='replace', action='write')
      log_unit = unit
      log_unit_error = unit
    end if

    call execute_forward(main_file=main_file, para_file=para_file, out_file=out_file, cwd=cwd)

    ! Complete Fortran output before returning control to Python.
    flush(log_unit)
    if (log_unit_error /= log_unit) flush(log_unit_error)
    if (len_trim(log_file) > 0) close(unit)
    log_unit = saved_unit
    log_unit_error = saved_error_unit
  end subroutine run_forward
end module runner
