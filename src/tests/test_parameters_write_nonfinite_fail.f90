program test_parameters_write_nonfinite_fail
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use mo_kind, only: dp
  use mo_main_config, only: parameters_t
  use mo_parameter_types, only: parameter_t

  implicit none

  type(parameters_t) :: parameters
  type(parameter_t) :: values(1)
  real(dp) :: invalid(1)

  values(1) = parameter_t(value=1.0_dp, min=0.0_dp, max=2.0_dp)
  invalid(1) = ieee_value(invalid(1), ieee_quiet_nan)
  call parameters%begin_configuration()
  call parameters%add_process("pet", values, [character(8) :: "factor"], group="pet_m1")
  call parameters%seal()
  call parameters%write_namelist("invalid_parameter_nonfinite.nml", invalid)

  error stop "non-finite parameter output was accepted"
end program test_parameters_write_nonfinite_fail
