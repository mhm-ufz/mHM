!> \file    mo_parameter_namelists.F90
!> \copydoc mo_parameter_namelists

!> \brief   Central input configuration for process parameters.
!> \details This type owns the generated parameter namelist instances exposed
!! through the exchange configuration. Process containers remain responsible
!! for selecting, loading, validating, and registering the blocks they use.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Jul 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_parameter_namelists
  use mo_logging
  use nml_helper, only: NML_OK
  use nml_interception_1, only: nml_interception_1_t
  use nml_snow_1, only: nml_snow_1_t
  use nml_soil_moisture_1, only: nml_soil_moisture_1_t
  use nml_soil_moisture_2, only: nml_soil_moisture_2_t
  use nml_soil_moisture_3, only: nml_soil_moisture_3_t
  use nml_soil_moisture_4, only: nml_soil_moisture_4_t
  use nml_direct_runoff_1, only: nml_direct_runoff_1_t
  use nml_pet_m2, only: nml_pet_m2_t
  use nml_pet_m1, only: nml_pet_m1_t
  use nml_pet_1, only: nml_pet_1_t
  use nml_pet_2, only: nml_pet_2_t
  use nml_pet_3, only: nml_pet_3_t
  use nml_interflow_1, only: nml_interflow_1_t
  use nml_percolation_1, only: nml_percolation_1_t
  use nml_baseflow_1, only: nml_baseflow_1_t
  use nml_neutrons_1, only: nml_neutrons_1_t
  use nml_neutrons_2, only: nml_neutrons_2_t
  use nml_routing_1, only: nml_routing_1_t
  use nml_routing_2, only: nml_routing_2_t
  use nml_routing_3, only: nml_routing_3_t
  use nml_river_temperature_1, only: nml_river_temperature_1_t

  implicit none
  private

  !> \class parameter_namelists_t
  !> \brief Generated parameter input blocks for one domain.
  type, public :: parameter_namelists_t
    type(nml_interception_1_t) :: interception_1
    type(nml_snow_1_t) :: snow_1
    type(nml_soil_moisture_1_t) :: soil_moisture_1
    type(nml_soil_moisture_2_t) :: soil_moisture_2
    type(nml_soil_moisture_3_t) :: soil_moisture_3
    type(nml_soil_moisture_4_t) :: soil_moisture_4
    type(nml_direct_runoff_1_t) :: direct_runoff_1
    type(nml_pet_m2_t) :: pet_m2
    type(nml_pet_m1_t) :: pet_m1
    type(nml_pet_1_t) :: pet_1
    type(nml_pet_2_t) :: pet_2
    type(nml_pet_3_t) :: pet_3
    type(nml_interflow_1_t) :: interflow_1
    type(nml_percolation_1_t) :: percolation_1
    type(nml_baseflow_1_t) :: baseflow_1
    type(nml_neutrons_1_t) :: neutrons_1
    type(nml_neutrons_2_t) :: neutrons_2
    type(nml_routing_1_t) :: routing_1
    type(nml_routing_2_t) :: routing_2
    type(nml_routing_3_t) :: routing_3
    type(nml_river_temperature_1_t) :: river_temperature_1
  contains
    procedure :: set_dims => parameter_namelists_set_dims
  end type parameter_namelists_t

contains

  !> \brief Set runtime dimensions for generated parameter namelists.
  subroutine parameter_namelists_set_dims(self, n_geo_units)
    class(parameter_namelists_t), intent(inout), target :: self
    integer, intent(in) :: n_geo_units
    character(1024) :: errmsg
    integer :: status

    status = self%baseflow_1%set_dims(n_geo_units=n_geo_units, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Error setting baseflow parameter dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine parameter_namelists_set_dims

end module mo_parameter_namelists
