program test_parameters_strict_bounds_fail
  use mo_kind, only: dp
  use mo_main_config, only: parameters_t
  use nml_helper, only: parameter_t

  implicit none

  type(parameters_t) :: parameters
  type(parameter_t) :: values(1)

  values(1) = parameter_t(value=1.0_dp, optimize=.true., lower_bound=1.0_dp, upper_bound=1.0_dp)
  call parameters%begin_configuration()
  call parameters%add_process("pet", values, [character(8) :: "factor"], group="pet_m1")
  call parameters%seal()
  error stop "optimized parameter with non-strict bounds was accepted"
end program test_parameters_strict_bounds_fail
