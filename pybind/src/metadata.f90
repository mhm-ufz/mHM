!> \file metadata.f90
!> \copydoc metadata

!> \brief Compiled version and forward-run defaults exposed through f2py.
!> \details Uses the same version and path constants as the native driver.
!> \authors Sebastian Mueller
!> \date Nov 2022
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup py_mhm
module metadata
  use mo_version, only: compiled_version => version
  implicit none
contains
  !> \brief Return the compiled Fortran model version.
  function version() result(value)
    character(256) :: value
    value = compiled_version
  end function version

  !> \brief Return the canonical default paths without duplicating them in Python.
  subroutine run_defaults(main_file, para_file, out_file, cwd)
    use mo_common_constants, only: default_main_file, default_para_file, default_out_file, default_cwd
    character(256), intent(out) :: main_file !< Default main namelist path.
    character(256), intent(out) :: para_file !< Default parameter namelist path.
    character(256), intent(out) :: out_file !< Default output namelist path.
    character(256), intent(out) :: cwd !< Default working directory.
    main_file = default_main_file
    para_file = default_para_file
    out_file = default_out_file
    cwd = default_cwd
  end subroutine run_defaults
end module metadata
