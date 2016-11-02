!> \file mrm_driver.f90
#ifndef MRM2MHM
program mrm_driver

  use mo_common_variables,              only: global_parameters, global_parameters_name, optimize
  use mo_finish,                        only: finish
  use mo_kind,                          only: dp, i4
  use mo_message,                       only: message
  use mo_mrm_eval,                      only: mrm_eval
  use mo_mrm_global_variables,          only: dirConfigOut
  use mo_mrm_init,                      only: mrm_init
  use mo_mrm_objective_function_runoff, only: single_objective_runoff
  use mo_mrm_write,                     only: mrm_write, mrm_write_optifile, mrm_write_optinamelist
  use mo_optimization,                  only: optimization
  use mo_string_utils,                  only: num2str
  use mo_timer,                         only: timers_init, timer_start, timer_stop, timer_get

  
  implicit none
  
  ! variables for optimization
  real(dp)             :: funcbest    ! best objective function achivied during optimization
  integer(i4)          :: itimer
  logical, allocatable :: maskpara(:) ! true  = parameter will be optimized     = parameter(i,4) = 1
  !                                   ! false = parameter will not be optimized = parameter(i,4) = 0
  ! --------------------------------------------------------------------------
  ! INITIALIZE
  ! --------------------------------------------------------------------------
  call mrm_init()
  ! Start timings
  print*, 'start timer'
  call timers_init
  if (optimize) then
     ! -----------------------------------------------------------------------
     ! OPTIMIZE
     ! -----------------------------------------------------------------------
     call optimization(single_objective_runoff, dirConfigOut, funcbest, maskpara)
     ! write a file with final objective function and the best parameter set
     call mrm_write_optifile(funcbest, global_parameters(:,3), global_parameters_name(:))
     ! write a file with final best parameter set in a namelist format
     call mrm_write_optinamelist(global_parameters, maskpara, global_parameters_name(:))
     deallocate(maskpara)
  else
     ! -----------------------------------------------------------------------
     ! FORWARD RUN
     ! -----------------------------------------------------------------------
     itimer = 1
     call timer_start(itimer)
     call message('    perform forward run of mRM')
     call mrm_eval(global_parameters(:,3))
     call timer_stop(itimer)
     call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  end if
  ! --------------------------------------------------------------------------
  ! WRITE OUTPUT
  ! --------------------------------------------------------------------------
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
