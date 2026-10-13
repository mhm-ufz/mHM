!> \brief Version-only wrapper retained for future Python bindings.
!> \authors Sebastian Mueller
!> \copyright Copyright 2005-today, the mHM Developers. LGPLv3+.
module metadata
  use mo_version, only: compiled_version => version
  implicit none
contains
  !> \brief Return the compiled Fortran model version.
  function version() result(value)
    character(256) :: value
    value = compiled_version
  end function version
end module metadata
