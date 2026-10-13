!> \file    mo_geology_classdefinition.f90
!> \brief   Definition for geology class lookup container.
!> \details Provides a light-weight container that stores information read from
!!          `geology_classdefinition.txt`. The container offers a `reset`
!!          procedure to clear allocated memory when the data needs to be
!!          re-used.
module mo_geology_classdefinition

  use mo_kind, only : i4

  implicit none

  private
  public :: geology_classdefinition_t

  !> \class   geology_classdefinition_t
  !> \brief   Holds geological class lookup data.
  type :: geology_classdefinition_t
    integer(i4) :: nGeo = 0_i4                         !< number of geological classes
    integer(i4), allocatable :: geo_param_index(:)      !< global GeoParam index per geological class
    integer(i4), allocatable :: geo_unit(:)            !< geological unit identifiers
    integer(i4), allocatable :: geo_karstic(:)         !< karstic flag per geological unit
  contains
    procedure, public :: reset => geology_classdefinition_reset
  end type geology_classdefinition_t

contains

  !> \brief Reset the geology class definition container.
  !> \details Deallocates all allocatable members and resets counters.
  subroutine geology_classdefinition_reset(self)
    class(geology_classdefinition_t), intent(inout) :: self

    if (allocated(self%geo_param_index)) deallocate(self%geo_param_index)
    if (allocated(self%geo_unit)) deallocate(self%geo_unit)
    if (allocated(self%geo_karstic)) deallocate(self%geo_karstic)
    self%nGeo = 0_i4
  end subroutine geology_classdefinition_reset

end module mo_geology_classdefinition
