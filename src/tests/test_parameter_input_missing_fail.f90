program test_parameter_input_missing_fail
  use mo_domain, only: domain_t

  implicit none

  type(domain_t), target :: domain

  call domain%create(main_file="test_nml/mrm_minimal.nml", cwd=".")
  call domain%configure(main_file="test_nml/mrm_minimal.nml")

end program test_parameter_input_missing_fail
