!> \file mo_optimization_types.f90

!> \authors Maren Kaluza
!> \date November 2019

MODULE mo_optimization_types
  use mo_kind, only : i4, dp

  ! Written Maren Kaluza, November 2019

  IMPLICIT NONE

  public :: optidata

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
  end type optidata_sim
  

END MODULE mo_optimization_types
