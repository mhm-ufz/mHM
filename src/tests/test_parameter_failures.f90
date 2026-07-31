program test_parameter_failures
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use, intrinsic :: iso_fortran_env, only: error_unit
  use mo_domain, only: domain_t
  use mo_kind, only: dp
  use mo_main_config, only: parameters_t
  use mo_parameter_types, only: parameter_t

  implicit none

  character(64) :: scenario

  if (command_argument_count() /= 1) then
    write(error_unit, '(a)') "usage: test_parameter_failures <scenario>"
    error stop 2
  end if
  call get_command_argument(1, scenario)

  select case (trim(scenario))
  case ("duplicate-name")
    call duplicate_name()
  case ("duplicate-group")
    call duplicate_group()
  case ("definition-bounds")
    call definition_bounds()
  case ("initial-bounds")
    call initial_bounds()
  case ("strict-bounds")
    call strict_bounds()
  case ("set-bounds")
    call set_bounds()
  case ("write-bounds")
    call write_bounds()
  case ("write-nonfinite")
    call write_nonfinite()
  case ("write-size")
    call write_size()
  case ("missing-input")
    call missing_input()
  case default
    write(error_unit, '(a)') "unknown parameter failure scenario: " // trim(scenario)
    error stop 2
  end select

  write(error_unit, '(a)') "parameter failure scenario unexpectedly succeeded: " // trim(scenario)

contains

  subroutine duplicate_name()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(2)

    definitions = parameter_t(value=1.0_dp, min=0.0_dp, max=2.0_dp)
    call parameters%begin_configuration()
    call parameters%add_process("pet", definitions, [character(8) :: "factor", "FACTOR"], group="pet_m1")
  end subroutine duplicate_name

  subroutine duplicate_group()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(1)

    definitions = parameter_t(value=1.0_dp, min=0.0_dp, max=2.0_dp)
    call parameters%begin_configuration()
    call parameters%add_process("pet", definitions, [character(8) :: "factor"], group="pet_m1")
    call parameters%add_process("other", definitions, [character(8) :: "factor"], group="PET_M1")
  end subroutine duplicate_group

  subroutine definition_bounds()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(1)

    definitions = parameter_t(value=1.0_dp, optimize=.true., min=2.0_dp, max=0.0_dp)
    call configure_and_seal(parameters, definitions)
  end subroutine definition_bounds

  subroutine initial_bounds()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(1)

    definitions = parameter_t(value=3.0_dp, min=0.0_dp, max=2.0_dp)
    call configure_and_seal(parameters, definitions)
  end subroutine initial_bounds

  subroutine strict_bounds()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(1)

    definitions = parameter_t(value=1.0_dp, optimize=.true., min=1.0_dp, max=1.0_dp)
    call configure_and_seal(parameters, definitions)
  end subroutine strict_bounds

  subroutine set_bounds()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(1)

    definitions = parameter_t(value=1.0_dp, min=0.0_dp, max=2.0_dp)
    call configure_and_seal(parameters, definitions)
    call parameters%set([3.0_dp])
  end subroutine set_bounds

  subroutine write_bounds()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(1)

    definitions = parameter_t(value=1.0_dp, min=0.0_dp, max=2.0_dp)
    call configure_and_seal(parameters, definitions)
    call parameters%write_namelist("invalid_parameter_bounds.nml", [3.0_dp])
  end subroutine write_bounds

  subroutine write_nonfinite()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(1)
    real(dp) :: invalid(1)

    definitions = parameter_t(value=1.0_dp, min=0.0_dp, max=2.0_dp)
    invalid(1) = ieee_value(invalid(1), ieee_quiet_nan)
    call configure_and_seal(parameters, definitions)
    call parameters%write_namelist("invalid_parameter_nonfinite.nml", invalid)
  end subroutine write_nonfinite

  subroutine write_size()
    type(parameters_t) :: parameters
    type(parameter_t) :: definitions(1)

    definitions = parameter_t(value=1.0_dp, min=0.0_dp, max=2.0_dp)
    call configure_and_seal(parameters, definitions)
    call parameters%write_namelist("invalid_parameter_size.nml", [1.0_dp, 1.5_dp])
  end subroutine write_size

  subroutine missing_input()
    type(domain_t), target :: domain

    call domain%create(main_file="test_nml/mrm_minimal.nml", cwd=".")
    call domain%configure(main_file="test_nml/mrm_minimal.nml")
  end subroutine missing_input

  subroutine configure_and_seal(parameters, definitions)
    type(parameters_t), intent(out) :: parameters
    type(parameter_t), intent(in) :: definitions(:)

    call parameters%begin_configuration()
    call parameters%add_process("pet", definitions, [character(8) :: "factor"], group="pet_m1")
    call parameters%seal()
  end subroutine configure_and_seal

end program test_parameter_failures
