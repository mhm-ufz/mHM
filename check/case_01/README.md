These are some instructions for recreating the restart files in case of back-breaking changes.

Create an `mhm.nml` to create the reference datasets.  
```
    write_restart = .true.
    read_restart = .false.
    dir_out(1) = 'output_save/b1_'
    mhm_file_restartout(1) = 'output_save/b1_mHM_restart_001.nc'
    mrm_file_restartout(1) = 'output_save/b1_mRM_restart_001.nc'
    dirconfigout = 'output_save/'
    eval_per(1)%dend = 31
    eval_per(1)%dstart = 2
    eval_per(1)%mend = 12
    eval_per(1)%mstart = 7
    eval_per(1)%yend = 1990
    eval_per(1)%ystart = 1990
    warming_days(1) = 182
```
Create an `mhm.nml` to create the restart files.  
```
    write_restart = .true.
    read_restart = .false.
    dir_out(1) = 'restart/b1_'
    mhm_file_restartout(1) = 'restart/b1_mHM_restart_001.nc'
    mrm_file_restartout(1) = 'restart/b1_mRM_restart_001.nc'
    dirconfigout = 'restart/'
    eval_per(1)%dend = 1
    eval_per(1)%dstart = 1
    eval_per(1)%mend = 7
    eval_per(1)%mstart = 7
    eval_per(1)%yend = 1990
    eval_per(1)%ystart = 1990
    warming_days(1) = 181
```
Use this `mhm.nml` for the check.  
```
    write_restart = .true.
    read_restart = .true.
    dir_out(1) = 'output_b1/b1_'
    mhm_file_restartout(1) = 'output_b1/b1_mHM_restart_001.nc'
    mrm_file_restartout(1) = 'output_b1/b1_mRM_restart_001.nc'
    dirconfigout = 'output_b1/'
    eval_per(1)%dend = 31
    eval_per(1)%dstart = 2
    eval_per(1)%mend = 12
    eval_per(1)%mstart = 7
    eval_per(1)%yend = 1990
    eval_per(1)%ystart = 1990
    warming_days(1) = 0
```