!> \file mrm_driver.f90
#ifndef mrm2mhm
program mrm_driver
  use mo_mrm_init, only: mrm_init
  use mo_mrm_eval, only: mrm_eval
  use mo_mrm_write, only: mrm_write
  use mo_finish, only: finish
  use mo_mrm_global_variables, only: mrm_global_parameters
  implicit none
  call mrm_init()
  call mrm_eval(mrm_global_parameters(:,3))
  call mrm_write()
  ! --------------------------------------------------------------------------
  ! FINISH UP
  ! --------------------------------------------------------------------------
  call finish('mRM','Finished!')
end program mrm_driver
#else
! dummy module such that this file is never empty for compilation
module dummy
  implicit none
end module dummy
#endif
