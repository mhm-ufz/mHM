!> \file    app/driver.F90
!> \copydoc driver

!> \brief mHM v6 command-line driver.
!> \details Parses command-line options, configures logging and message visibility,
!! and delegates forward execution to mo_mhm_runner.
!> \authors Sebastian Mueller
!> \date    Jan 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
#include "logging.h"
program driver
  use mo_logging
  use mo_cli, only: cli_parser
  use mo_message, only: use_log
  use mo_cli_support, only: set_verbosity_level
  use mo_version, only: version
  use mo_kind, only: i4
  use mo_common_constants, only: default_main_file, default_para_file, default_out_file, default_cwd
  use mo_mhm_runner, only: run_forward
  implicit none

  type(cli_parser) :: parser
  integer :: unit

  parser = cli_parser(prog="mhm", description="The mesoscale hydrological model - mHM v6", &
    add_logger_options=.true., add_version_option=.true., version=version)

  call parser%add_option(name='log-file', has_value=.true., required=.false., &
    value_name="path", help='Write to log file.')

  call parser%add_option(name="nml", s_name="n", has_value=.true., &
    value_name="path", default=default_main_file, help="The mHM configuration namelist.")

  call parser%add_option(name="parameters", s_name="p", has_value=.true., &
    value_name="path", default=default_para_file, help="The mHM parameter namelist.")

  call parser%add_option(name="output", s_name="o", has_value=.true., &
    value_name="path", default=default_out_file, help="The mHM output namelist.")

  call parser%add_option(name="cwd", blank=.true., default=default_cwd, help="The desired working directory (optional).")
  call parser%parse()

  use_log = .true. ! Enable message logging after parsing so help/version output stays clean.
  if (parser%option_was_read('log-file')) then
    open(newunit=unit, file=parser%option_value('log-file'), status='replace', action='write')
    log_unit = unit
    log_unit_error = unit
  end if
  call set_verbosity_level(3_i4 - parser%option_read_count("quiet"))

  call run_forward(main_file=parser%option_value("nml"), para_file=parser%option_value("parameters"), &
    out_file=parser%option_value("output"), cwd=parser%option_value("cwd"))
end program driver
