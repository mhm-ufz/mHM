!> \file mo_optimization_types.f90

!> \authors Maren Kaluza
!> \date November 2019

MODULE mo_optimization_types
  use mo_kind, only : i4, dp

  ! Written Maren Kaluza, November 2019

  IMPLICIT NONE

  public :: optidata, optidata_sim

  private

  ! optional data, such as sm, neutrons, et, tws
  ! data type for observed data, providing metadata
  ! for simulated data
  ! dim1 = number grid cells L1
  ! dim2 = number of meteorological time steps
  type optidata
    real(dp), dimension(:, :), allocatable    :: dataObs ! observed data
    logical, dimension(:, :), allocatable     :: maskObs ! mask of observed data
    character(256)                            :: dir ! directory where to read opti data
    integer(i4)                               :: timeStepInput ! time step of optional data
    integer(i4)                               :: writeOutCounter ! the current timestep
                                                                 ! the simulated opti data is written to
  end type optidata

  ! type for simulated optional data
  type optidata_sim
    real(dp), dimension(:, :), allocatable    :: dataSim
    integer(i4)                               :: writeOutCounter ! the current timestep
                                                                 ! the simulated opti data is written to
    contains
    procedure :: init => optidata_sim_init
    procedure :: destroy => optidata_sim_destroy
  end type optidata_sim

  contains

  subroutine optidata_sim_init(this, optidataObs)
    class(optidata_sim), intent(inout) :: this
    type(optidata),      intent(in)    :: optidataObs

    allocate(this%dataSim(size(optidataObs%dataObs, dim = 1), size(optidataObs%dataObs, dim = 2)))
    this%dataSim(:, :) = 0.0_dp ! has to be intialized with zero because later summation
    this%writeOutCounter = 1
  end subroutine optidata_sim_init

  subroutine optidata_sim_destroy(this)
    class(optidata_sim), intent(inout) :: this

    deallocate(this%dataSim)
  end subroutine optidata_sim_destroy

END MODULE mo_optimization_types
