!> \file test_forward_run.f90
!> \brief Check default arguments, path overrides, and working-directory isolation of forward execution.
!> \date Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
program test_forward_run
  use mo_mhm_runner, only: run_forward
  use mo_kind, only: dp
  use mo_os, only: path_cwd
  use mo_netcdf, only: NcDataset, NcVariable
  implicit none

  character(:), allocatable :: original_cwd, absent
  real(dp), allocatable :: reference(:, :, :)

  original_cwd = path_cwd()
  call run_forward()
  call read_discharge("mrm_output.nc", reference)

  ! An alternate run root must resolve all default namelists and input/output paths.
  call run_forward(cwd="alternate")
  call check_output("alternate/mrm_output.nc")

  ! Explicit names and cwd must match default-name execution in the same root.
  call run_forward(main_file="test_nml/mrm_minimal.nml", para_file="test_nml/mhm_parameter_v6_routing.nml", &
    out_file="mhm-output-template.nml", cwd="alternate")
  call check_output("alternate/mrm_output.nc")

  ! Forwarding unallocated optional actuals must select defaults on a repeated invocation.
  call run_forward(main_file=absent, para_file=absent, out_file=absent, cwd=absent)
  call check_output("mrm_output.nc")

contains

  !> \brief Read the routed discharge produced by a forward run.
  subroutine read_discharge(file, discharge)
    character(*), intent(in) :: file !< Gridded output path.
    real(dp), allocatable, intent(out) :: discharge(:, :, :) !< Discharge records on the routing grid.
    type(NcDataset) :: nc
    type(NcVariable) :: var

    nc = NcDataset(file, "r")
    var = nc%getVariable("discharge")
    call var%getData(discharge)
    call nc%close()
  end subroutine read_discharge

  !> \brief Compare output records and verify the process working directory was preserved.
  subroutine check_output(file)
    character(*), intent(in) :: file !< Forward-run output to compare against default execution.
    real(dp), allocatable :: actual(:, :, :)

    if (path_cwd() /= original_cwd) error stop "forward run changed the process working directory"
    call read_discharge(file, actual)
    if (any(shape(actual) /= shape(reference))) error stop "forward run changed the discharge output shape"
    if (any(actual /= reference)) error stop "forward run path/default handling changed discharge values"
  end subroutine check_output

end program test_forward_run
