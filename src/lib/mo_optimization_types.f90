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
    character(256)                            :: varname
  end type optidata

  ! type for simulated optional data
  type optidata_sim
    real(dp), dimension(:, :), allocatable    :: dataSim
    integer(i4)                               :: writeOutCounter ! the current timestep
                                                                 ! the simulated opti data is written to
    integer(i4)                               :: averageCounter  ! set to 0 on average, incremented on add

    contains
    procedure :: init => optidata_sim_init
    procedure :: destroy => optidata_sim_destroy
    procedure :: increment_counter => optidata_sim_increment_counter
    procedure :: add => optidata_sim_add
    procedure :: average => optidata_sim_average
    procedure :: average_per_timestep => optidata_sim_average_per_timestep
    procedure :: average_add => optidata_sim_average_add
  end type optidata_sim

  contains

  subroutine optidata_sim_init(this, optidataObs)
    class(optidata_sim), intent(inout) :: this
    type(optidata),      intent(in)    :: optidataObs

    allocate(this%dataSim(size(optidataObs%dataObs, dim = 1), size(optidataObs%dataObs, dim = 2)))
    this%dataSim(:, :) = 0.0_dp ! has to be intialized with zero because later summation
    this%writeOutCounter = 1
    this%averageCounter = 0
  end subroutine optidata_sim_init

  subroutine optidata_sim_destroy(this)
    class(optidata_sim), intent(inout) :: this

    deallocate(this%dataSim)
  end subroutine optidata_sim_destroy

  subroutine optidata_sim_increment_counter(this, timeStepInput, is_new_day, is_new_month, is_new_year)
    class(optidata_sim), intent(inout) :: this
    integer(i4),         intent(in)    :: timeStepInput
    logical,             intent(in)    :: is_new_day
    logical,             intent(in)    :: is_new_month
    logical,             intent(in)    :: is_new_year

    select case(timeStepInput)
    case(-1) ! daily
      if (is_new_day)   then
        this%writeOutCounter = this%writeOutCounter + 1
      end if
    case(-2) ! monthly
      if (is_new_month) then
        this%writeOutCounter = this%writeOutCounter + 1
      end if
    case(-3) ! yearly
      if (is_new_year)  then
        this%writeOutCounter = this%writeOutCounter + 1
      end if
    end select

  end subroutine optidata_sim_increment_counter

  subroutine optidata_sim_add(this, data_sim)
    class(optidata_sim),    intent(inout) :: this
    real(dp), dimension(:), intent(in)    :: data_sim

    this%dataSim(:, this%writeOutCounter) = &
            this%dataSim(:, this%writeOutCounter) + data_sim(:)
  end subroutine optidata_sim_add

  subroutine optidata_sim_average(this)
    class(optidata_sim), intent(inout) :: this

    this%dataSim(:, this%writeOutCounter) = &
            this%dataSim(:, this%writeOutCounter) / real(this%averageCounter, dp)
    this%writeOutCounter = this%writeOutCounter + 1
    this%averageCounter = 0
  end subroutine optidata_sim_average

  subroutine optidata_sim_average_per_timestep(this, timeStepInput, is_new_day, is_new_month, is_new_year)
    class(optidata_sim), intent(inout) :: this
    integer(i4),         intent(in)    :: timeStepInput
    logical,             intent(in)    :: is_new_day
    logical,             intent(in)    :: is_new_month
    logical,             intent(in)    :: is_new_year

    select case(timeStepInput)
    case(-1) ! daily
      if (is_new_day)   then
        call this%average()
      end if
    case(-2) ! monthly
      if (is_new_month) then
        call this%average()
      end if
    case(-3) ! yearly
      if (is_new_year)  then
        call this%average()
      end if
    end select
  end subroutine optidata_sim_average_per_timestep

  subroutine optidata_sim_average_add(this, data_sim)
    class(optidata_sim),    intent(inout) :: this
    real(dp), dimension(:), intent(in)    :: data_sim

    call this%add(data_sim(:))
    this%averageCounter = this%averageCounter + 1
  end subroutine optidata_sim_average_add

END MODULE mo_optimization_types
