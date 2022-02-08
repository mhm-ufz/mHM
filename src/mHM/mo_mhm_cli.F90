!> \file    mo_mhm_cli.f90
!> \brief   Module to parse command line arguments of mHM.
!> \details \copydetails mo_mhm_cli

!> \brief   Module to parse command line arguments of mHM.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Oct 2021
!> \details A simple parser for command line arguments for mHM.
!!          You can pass the path to mhm.nml, mhm_parameters.nml, mhm_outputs.nml and mrm_outputs.nml
!!
!!          You can also pass the CWD as plain last argument and get a help or version text.
!!
!!          \code{.sh}
!!          The mesoscale hydrological model - mHM
!!
!!          Usage: mhm [options] <cwd>
!!
!!          Options:
!!          <cwd>
!!              Description: The desired working directory (optional).
!!
!!          --help / -h
!!              Description: Print this help message.
!!
!!          --version / -v
!!              Description: Print the version of the program.
!!
!!          --nml / -n <path>
!!              Description: The mHM configuration namelist.
!!              Default: mhm.nml
!!
!!          --parameter / -p <path>
!!              Description: The mHM parameter namelist.
!!              Default: mhm_parameter.nml
!!
!!          --mhm_output / -o <path>
!!              Description: The mHM output namelist.
!!              Default: mhm_output.nml
!!
!!          --mrm_output / -r <path>
!!              Description: The mRM output namelist.
!!              Default: mrm_output.nml
!!          \endcode
module mo_mhm_cli

#ifdef NAG
  USE f90_unix_dir, ONLY : CHDIR
#endif

  implicit none

  private

  public :: parse_command_line

contains

  !> \brief parse the given command line arguments.
  subroutine parse_command_line()
    use mo_cli, only: cli_parser
    use mo_file, only: version, file_namelist_mhm, file_namelist_mhm_param, file_defOutput
    use mo_mrm_file, only: mrm_file_defOutput => file_defOutput

    implicit none

    type(cli_parser) :: parser

    parser = cli_parser( &
      description="The mesoscale hydrological model - mHM", &
      add_version_option=.true., &
      version=version)

    call parser%add_option( &
      name="cwd", &
      blank=.true., &
      help="The desired working directory (optional).")

    call parser%add_option( &
      name="nml", &
      s_name="n", &
      has_value=.true., &
      value_name="path", &
      default="mhm.nml", &
      help="The mHM configuration namelist.")

    call parser%add_option( &
      name="parameter", &
      s_name="p", &
      has_value=.true., &
      value_name="path", &
      default="mhm_parameter.nml", &
      help="The mHM parameter namelist.")

    call parser%add_option( &
      name="mhm_output", &
      s_name="o", &
      has_value=.true., &
      value_name="path", &
      default="mhm_outputs.nml", &
      help="The mHM output namelist.")

    call parser%add_option( &
      name="mrm_output", &
      s_name="r", &
      has_value=.true., &
      value_name="path", &
      default="mrm_outputs.nml", &
      help="The mRM output namelist.")

    ! parse given command line arguments
    call parser%parse()

    ! change working directory first
    if (parser%option_was_read("cwd")) call chdir(parser%option_value("cwd"))

    ! set nml file paths
    file_namelist_mhm = parser%option_value("nml")
    file_namelist_mhm_param = parser%option_value("parameter")
    file_defOutput = parser%option_value("mhm_output")
    mrm_file_defOutput = parser%option_value("mrm_output")

  end subroutine parse_command_line

end module mo_mhm_cli
