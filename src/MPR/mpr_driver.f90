!> \file mpr_driver.f90
#ifdef MPR_STANDALONE
program mpr_driver

  use mo_finish, only : finish
  use mo_mpr_eval, only : mpr_eval
  use mo_read_wrapper, only : read_data
  use mo_mpr_read_config, only : mpr_read_config
  USE mo_common_read_config, ONLY : common_read_config                    ! Read main configuration files
  use mo_common_variables, only : dirRestartOut, write_restart
  use mo_timer, only : timers_init
  use mo_mpr_startup, only : mpr_initialize
  use mo_mpr_restart, only : write_mpr_restart_files
  use mo_mpr_file, only : &
          file_namelist_mpr_param, unamelist_mpr_param, &      ! filename of namelist: mhm model parameter
          file_namelist_mpr, unamelist_mpr ! file containing main configurations
  use mo_mpr_global_variables, only : c2TSTu

  implicit none

  ! --------------------------------------------------------------------------
  ! INITIALIZE
  ! --------------------------------------------------------------------------
  call common_read_config(file_namelist_mpr, unamelist_mpr)
  call mpr_read_config(file_namelist_mpr, unamelist_mpr, file_namelist_mpr_param, unamelist_mpr_param)

  ! Start timings
  call timers_init

  call read_data()

  ! TODO: this might become part of the namelist in MPR-STANDALONE?!
  c2TSTu = 1.0_dp / 24.0_dp
  call mpr_initialize()

  ! -----------------------------------------------------------------------
  ! EXECUTION
  ! -----------------------------------------------------------------------
  call mpr_eval()

  ! --------------------------------------------------------------------------
  ! WRITE OUTPUT
  ! --------------------------------------------------------------------------
  if (write_restart) then
    call write_mpr_restart_files(dirRestartOut)
  end if
  ! --------------------------------------------------------------------------
  ! FINISH UP
  ! --------------------------------------------------------------------------
  call finish('MPR', 'Finished!')
end program mpr_driver
#else
! dummy module such that this file is never empty for compilation
module dummy_mpr
  implicit none
end module dummy_mpr
#endif
