program test_parameters_write_size_fail
  use mo_kind, only: dp
  use mo_main_config, only: parameters_t
  use mo_parameter_types, only: parameter_t

  implicit none

  type(parameters_t) :: parameters
  type(parameter_t) :: values(1)

  values(1) = parameter_t(value=1.0_dp, optimize=.false., min=0.0_dp, max=2.0_dp)
  call parameters%begin_configuration()
  call parameters%add_process("pet", values, [character(8) :: "factor"], group="pet_m1")
  call parameters%seal()
  call parameters%write_namelist("invalid_parameter_size.nml", [1.0_dp, 1.5_dp])

  error stop "wrong-sized parameter output was accepted"
end program test_parameters_write_size_fail
