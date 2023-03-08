!> \file mo_common_datetime_type.f90
!> \brief \copybrief mo_common_datetime_type
!> \details \copydetails mo_common_datetime_type

!> \brief type for date time information with an increment subroutine
!> \details Contains a current day, month, year, hour matching newTime, aswell as
!! previous day, month, year. Theses all get updated on increment
!! also contains nTimestep, and tIndex_out for writing
!! finally, contains iLAI and yId that are time dependent and updating routines
!! for these, and a function returning a boolean for writeout, dependent on the
!! timestep_model_input
!> \author Maren Kaluza
!> \date March 2019
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
MODULE mo_common_datetime_type
  use mo_kind, only : i4, dp

  ! Written Maren Kaluza, March 2019

  IMPLICIT NONE

  public :: datetimeinfo

  private

  type datetimeinfo
    !> number of timesteps in simulation period
    integer(i4)          :: nTimeSteps
    !> starts with simPer(iDomain)%julStart, then increments with
    !> julday(...)
    real(dp)             :: newTime

    !> current day
    integer(i4)          :: day
    !> current month
    integer(i4)          :: month
    !> current year
    integer(i4)          :: year
    ! hour is only local, only used for calculating other output
    integer(i4)          :: hour

    ! prev_day, prev_month, prev_year also are only used to store data to
    ! calculate is_new_day/month/year. These are only used locally
    integer(i4), private :: prev_day
    integer(i4), private :: prev_month
    integer(i4), private :: prev_year

    ! flags for stepping into new period
    logical              :: is_new_day
    logical              :: is_new_month
    logical              :: is_new_year

    integer(i4)          :: iLAI
    integer(i4)          :: yId

    ! for writing netcdf file
    integer(i4)          :: tIndex_out

    contains
    procedure :: init => datetimeinfo_init
    procedure :: increment => datetimeinfo_increment
    procedure :: update_LAI_timestep => datetimeinfo_update_LAI_timestep
    ! function, returns boolean dependent on tIndex_out and is_new_{period}
    procedure :: writeout => datetimeinfo_writeout
  end type datetimeinfo

  contains

  subroutine datetimeinfo_init(this, iDomain)
    use mo_common_mHM_mRM_variables, only : LCyearId, simPer, nTstepDay
    use mo_julian, only : caldat
    class(datetimeinfo), intent(inout) :: this
    integer(i4),     intent(in)    :: iDomain

    ! calculate NtimeSteps for this basin
    this%nTimeSteps = (simPer(iDomain)%julEnd - simPer(iDomain)%julStart + 1) * nTstepDay

    ! reinitialize time counter for LCover and MPR
    ! -0.5 is due to the fact that dec2date routine
    !   changes the day at 12:00 in NOON
    ! Whereas mHM needs day change at 00:00 h
    ! initialize the julian day as real
    this%newTime = real(simPer(iDomain)%julStart, dp)
    ! initialize the date
    call caldat(int(this%newTime), yy = this%year, mm = this%month, dd = this%day)
    ! initialize flags for period changes, they are true for first time step
    this%is_new_day   = .true.
    this%is_new_month = .true.
    this%is_new_year  = .true.

    ! initialize arrays and counters
    this%yId  = LCyearId(this%year, iDomain)
    this%hour = 0
    this%iLAI = 0

    ! this has no relevance yet. it is only so the variables are initialized
    this%prev_day   = this%day
    this%prev_month = this%month
    this%prev_year  = this%year

    this%tIndex_out = 0 ! tt if write out of warming period
  end subroutine datetimeinfo_init

  subroutine datetimeinfo_increment(this)
    use mo_julian, only : caldat, julday
    use mo_common_mHM_mRM_variables, only : timeStep
    class(datetimeinfo), intent(inout) :: this

    ! prepare the date and time information for next iteration step...
    ! set the current year as previous
    this%prev_day = this%day
    this%prev_month = this%month
    this%prev_year = this%year
    ! set the flags to false
    this%is_new_day   = .false.
    this%is_new_month = .false.
    this%is_new_year  = .false.

    ! increment of timestep
    this%hour = this%hour + timestep
    this%newTime = julday(this%day, this%month, this%year) + real(this%hour, dp) / 24._dp
    ! get correct hour for current day
    this%hour = mod(this%hour, 24)

    ! calculate new year, month and day
    call caldat(int(this%newTime), yy = this%year, mm = this%month, dd = this%day)
    ! update the flags
    if (this%prev_day   /= this%day) this%is_new_day = .true.
    if (this%prev_month /= this%month) this%is_new_month = .true.
    if (this%prev_year  /= this%year) this%is_new_year = .true.
  end subroutine datetimeinfo_increment

  subroutine datetimeinfo_update_LAI_timestep(this)
    use mo_mpr_global_variables, only : timeStep_LAI_input
    class(datetimeinfo), intent(inout) :: this

    select case (timeStep_LAI_input)
    case(0 : 1) ! long term mean monthly gridded fields or LUT-based values
      this%iLAI = this%month
    case(-1) ! daily timestep
      if (this%is_new_day) then
        this%iLAI = this%iLAI + 1
      end if
    case(-2) ! monthly timestep
      if (this%is_new_month) then
        this%iLAI = this%iLAI + 1
      end if
    case(-3) ! yearly timestep
      if (this%is_new_year) then
        this%iLAI = this%iLAI + 1
      end if
    end select
  end subroutine datetimeinfo_update_LAI_timestep

  function datetimeinfo_writeout(this, timeStep_model_outputs, tt) result(writeout)
    class(datetimeinfo), intent(in)    ::  this
    integer(i4),         intent(in)    ::  timeStep_model_outputs
    integer(i4),         intent(in)    ::  tt

    logical                            :: writeout

    writeout = .false.
    if (timeStep_model_outputs > 0) then
      if ((mod(this%tIndex_out, timeStep_model_outputs) == 0) &
                  .or. (tt == this%nTimeSteps)) writeout = .true.
    else
      select case(timeStep_model_outputs)
      case(0) ! only at last time step
        if (tt == this%nTimeSteps) writeout = .true.
      case(-1) ! daily
        if (((this%tIndex_out > 0) .and. this%is_new_day) .or. &
                  (tt == this%nTimeSteps))   writeout = .true.
      case(-2) ! monthly
        if (((this%tIndex_out > 0) .and. this%is_new_month) .or.&
                  (tt == this%nTimeSteps)) writeout = .true.
      case(-3) ! yearly
        if (((this%tIndex_out > 0) .and. this%is_new_year) .or. &
                  (tt == this%nTimeSteps))  writeout = .true.
      case default ! no output at all

      end select
    end if

  end function datetimeinfo_writeout

END MODULE mo_common_datetime_type
