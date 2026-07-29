program test_parameters_group_duplicate_fail
  use mo_kind, only: dp
  use mo_main_config, only: parameters_t
  use nml_helper, only: parameter_t

  implicit none

  type(parameters_t) :: parameters
  type(parameter_t) :: values(1)

  values(1) = parameter_t(value=1.0_dp, optimize=.false., lower_bound=0.0_dp, upper_bound=2.0_dp)
  call parameters%begin_configuration()
  call parameters%add_process("pet", values, [character(8) :: "factor"], group="pet_m1")
  call parameters%add_process("other", values, [character(8) :: "factor"], group="PET_M1")

  error stop "duplicate namelist group was accepted"
end program test_parameters_group_duplicate_fail
