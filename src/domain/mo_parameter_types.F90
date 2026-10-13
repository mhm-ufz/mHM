!> \file    mo_parameter_types.F90
!> \brief   \copybrief mo_parameter_types
!> \details \copydetails mo_parameter_types

!> \brief   Application-owned parameter configuration types.
!> \details The component order of parameter_t is part of the nml-tools
!! buffer-input contract and must match nml-schemas/parameter.yml.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jul 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
module mo_parameter_types
  use mo_kind, only: dp

  implicit none
  private

  !> \class parameter_t
  !> \brief Model parameter with optional calibration metadata.
  !> \details Component order is fixed as value, optimize, min, max.
  type, public :: parameter_t
    real(dp) :: value !< Parameter value
    logical :: optimize = .false. !< Include parameter in calibration
    real(dp) :: min !< Minimum calibration value
    real(dp) :: max !< Maximum calibration value
  end type parameter_t

end module mo_parameter_types
