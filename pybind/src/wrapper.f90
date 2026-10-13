!> \file wrapper.f90
!> \brief \copybrief metadata
!> \details \copydetails metadata

!> \brief Version-only wrapper retained for future Python bindings.
!> \details Exposes the shared compiled Fortran version through f2py.
!> \authors Sebastian Mueller
!> \date Nov 2022
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
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
