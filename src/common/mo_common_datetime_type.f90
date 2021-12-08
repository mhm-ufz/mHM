!> \file mo_common_datetime_type.f90

!< author: Maren Kaluza
!< date: March 2019
!< summary: types for datetime information, contains:
!<   - type: datetimeinfo
!<      - contains a current day, month, year, hour matching newTime, as well as previous day, month, year.
!<      - theses all get updated on increment routine
!<      - nTimestep, and tIndex_out for writing
!<
!< for these, and a function returning a boolean for writeout, dependent on the
!< timestep_model_input

MODULE mo_common_datetime_type
  use mo_kind, only : i4, dp, i8
  use mo_message, only: error_message, message
  use mo_string_utils, only: num2str

  ! Written Maren Kaluza, March 2019

  IMPLICIT NONE

  private

  public :: datetimeinfo, period, get_land_cover_period_indices

  !> [h] simulation time step (= TS) in [h]
  integer(i4), public :: timeStep
  !> Number of time intervals per day
  integer(i4), public :: nTstepDay

  type aggregateperiod
    !< multiple consecutive periods defining the coordinate information for mHM parameters
    !> number of periods/ids
    integer(i4) :: nIds
    !> name of periods (e.g. land cover or LAI)
    character(256) :: name
    !> currently active index (iterator incremented by increment procedure)
    integer(i4) :: i
    !> mapping of ids to each simPer year (TODO: merge into seondsSince/refJulDate structure at some point)
    integer(i4), dimension(:), allocatable :: yearIds
    !> all boundaries of input data (nIds + 1) in case of isRegular == .false. and timeStep yearly
    integer(i4), dimension(:), allocatable :: boundaries
    !> seconds since refJulDate (size nIds)
    integer(i8), dimension(:), allocatable :: secondsSince
    !> reference Julian date for time information
    integer(i4) :: refJulDate
    !> time step, by which aggregatePeriod is incremented
    integer(i4) :: timeStep
    !> flag showing, whether averaged values are used (day of year, month of year)
    logical :: isAveraged
    !> flag showing, whether timeStep is equal of all periods
    logical :: isRegular

  contains

    procedure, public :: init => init_aggregateperiod
    procedure :: increment => increment_aggregateperiod
    procedure :: get_values => get_values_aggregateperiod
    procedure :: get_unit => get_unit_aggregateperiod

  end type aggregateperiod

  !> periods for each domain specifying the lai time information
  type(aggregateperiod), dimension(:), allocatable, public :: laiPeriods
  !> periods for each domain specifying the land cover time information
  type(aggregateperiod), dimension(:), allocatable, public :: landCoverPeriods
  !> number of maximum lai periods defined globally through mhm.nml
  integer(i4), public :: nlaiPeriods
  !> number of maximum land cover periods defined globally through mhm.nml
  integer(i4), public :: nLandCoverPeriods

  type datetimeinfo
    !< type for describing a timestamp with additional information needed to handle mHM time loop
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
    integer(i4), private :: hour

    ! prev_day, prev_month, prev_year also are only used to store data to
    ! calculate is_new_day/month/year. These are only used locally
    integer(i4), private :: prev_day
    integer(i4), private :: prev_month
    integer(i4), private :: prev_year

    ! flags for stepping into new period
    logical              :: is_new_day
    logical              :: is_new_month
    logical              :: is_new_year

    ! for writing netcdf file
    integer(i4)          :: tIndex_out

    contains
    procedure :: init => datetimeinfo_init
    procedure :: increment => datetimeinfo_increment
    ! function, returns boolean dependent on tIndex_out and is_new_{period}
    procedure :: writeout => datetimeinfo_writeout
  end type datetimeinfo

  type period
    !< type for describing a period (time interval) defined by start and end date
    integer(i4) :: dStart      ! first day
    integer(i4) :: mStart      ! first month
    integer(i4) :: yStart      ! first year
    integer(i4) :: dEnd        ! last  day
    integer(i4) :: mEnd        ! last  month
    integer(i4) :: yEnd        ! last  year
    integer(i4) :: julStart    ! first julian day
    integer(i4) :: julEnd      ! last  julian day
    integer(i4) :: nObs        ! total number of observations

    contains

    procedure :: copy => period_copy_period_data
    procedure :: init => init_period
  end type period

  type(period), dimension(:), allocatable, public :: simPer      ! warmPer + evalPer

contains

  subroutine datetimeinfo_init(this, iDomain)
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

    ! initialize counters
    this%hour = -timestep

    ! this has no relevance yet. it is only so the variables are initialized
    this%prev_day   = this%day
    this%prev_month = this%month
    this%prev_year  = this%year

    this%tIndex_out = 1 ! tt if write out of warming period
  end subroutine datetimeinfo_init

  subroutine datetimeinfo_increment(this)
    use mo_julian, only : caldat, julday
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
    this%hour = mod(this%hour + timestep, 24)
    this%newTime = julday(this%day, this%month, this%year)&
            + real(this%hour + timestep, dp) / 24._dp
    ! calculate new year, month and day
    call caldat(int(this%newTime), yy = this%year, mm = this%month, dd = this%day)
    ! update the flags
    if (this%prev_day   /= this%day) this%is_new_day = .true.
    if (this%prev_month /= this%month) this%is_new_month = .true.
    if (this%prev_year  /= this%year) this%is_new_year = .true.
  end subroutine datetimeinfo_increment

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

  function get_doy(yy, mm, dd) result(doy)
    !< helper function returning the day of year for a given calendar day
    !> input day, month and year integer values
    ! TODO: think about moving to mo_julian
    integer(i4), intent(in)    ::  yy, mm, dd

    integer(i4), dimension(12), parameter :: months = [31,28,31,30,31,30,31,31,30,31,30,31]
    integer(i4) :: doy
    integer(i4) :: leapDay

    ! leap year
    leapDay = 0_i4
    if (is_leap_year(yy)) leapDay = 1_i4
    ! calculate doy as cumulative sum
    doy = sum(months(1:mm-1)) + dd + leapDay

  end function get_doy

  function is_leap_year(year) result(isLeapYear)
    !< \brief checks if a year is a leap year
    !< \details checks if a year is a leap year
    ! TODO: think about moving to mo_julian
    integer(i4), intent(in) :: year
    logical :: isLeapYear

    if (mod(year, 400_i4) == 0 ) then
      isLeapYear = .true.
    elseif (mod(year, 100_i4) == 0 ) then
      isLeapYear = .false.
    elseif (mod(year, 4_i4) == 0 ) then
      isLeapYear = .true.
    else
      isLeapYear = .false.
    end if
  end function is_leap_year

  subroutine period_copy_period_data(this, toPeriod)
    class(period), intent(inout) :: this
    class(period), intent(inout) :: toPeriod

    call toPeriod%init(this%dStart, this%mStart, this%yStart, this%dEnd, this%mEnd, this%yEnd)

  end subroutine period_copy_period_data

  subroutine init_period(this, dStart, mStart, yStart, dEnd, mEnd, yEnd)
    class(period), intent(inout) :: this
    integer(i4), intent(in) :: dStart, mStart, yStart, dEnd, mEnd, yEnd

    this%dStart = dStart   ! first day
    this%mStart = mStart    ! first month
    this%yStart = yStart    ! first year
    this%dEnd = dEnd      ! last  day
    this%mEnd = mEnd      ! last  month
    this%yEnd = yEnd      ! last  year
    this%julStart = 0 ! first julian day
    this%julEnd   = 0 ! last  julian day
    this%nObs     = 0 ! total number of observations

  end subroutine init_period

  subroutine init_aggregateperiod(this, n, nMax, name, simPerArg, units, periodValues, &
          keepUnneededPeriods, selectIndices)
    !< initializes an aggregateperiod instance

    class(aggregateperiod), intent(inout) :: this
    !> size of aggregatePeriod values and maximum size allowed (error is raised here, if violated)
    integer(i4), intent(in) :: n, nMax
    !> name of the aggregatePeriod (currently used for either land cover or LAI)
    character(*), intent(in) :: name
    !> reference simulation period (used for validity checking and selection)
    class(period), intent(in) :: simPerArg
    !> units of the values to describe aggregateperiod
    character(*), intent(in) :: units
    !> values for the aggregatePeriod (e.g. read from netcdf file)
    integer(i4), dimension(:), intent(in), optional :: periodValues
    !> flag whether to select periodValues so values not in simPerArg are discarded
    logical, intent(in), optional :: keepUnneededPeriods
    !> integer vector containing indices of periodValues to select (depends on keepUnneededPeriods)
    integer(i4), dimension(:), allocatable, intent(out) :: selectIndices

    integer(i4) :: iSel
    type(period) :: inPeriod

    this%nIds = n
    this%name = name
    ! sanity check
    if (n > nMax) then
      call error_message('There is a maximum of ', trim(num2str(nMax)), &
              trim(name), ' periods allowed. Passed ', trim(num2str(n)), '.')
    end if

    this%isRegular = .true.
    ! decide type of aggregateperiod based on units information
    select case(units)
    case('day_of_year', 'day of year')
      this%isAveraged = .true.
      this%timeStep = -1_i4
      ! keepUnneededPeriods is ignored (simPerArg is assumed to cover whole calender year), range 1:366 is assumed
      ! TODO: check this
      selectIndices = [(iSel, iSel=1, n)]
      ! select the doy for the simulation period start as initial counter
      this%i = get_doy(simPerArg%yStart, simPerArg%mStart, simPerArg%dStart)
    case('month_of_year', 'month of year')
      this%isAveraged = .true.
      this%timeStep = -2_i4
      ! keepUnneededPeriods is ignored (simPerArg is assumed to cover whole calender year), range 1:12 is assumed
      selectIndices = [(iSel, iSel=1, n)]
      ! select the moy for the simulation period start as initial counter
      this%i = simPerArg%mStart
    ! TODO: switch to a CF-compatible representation and naming
    case('years', 'year')
      ! this is custom land cover scene or land cover period information used in mHM
      if (.not. present(keepUnneededPeriods)) then
        call error_message('Cannot init aggregate period by units "', trim(units), &
                '", without flag keepUnneededPeriods.')
      end if
      this%isAveraged = .false.
      this%timeStep = -3_i4
      this%isRegular = .false.
      ! periodValues are boundaries here and thus are of size n+1 !!!
      call get_land_cover_period_indices(simPerArg, periodValues, keepUnneededPeriods, &
              yearIds=this%yearIds, selectIndices=selectIndices)
      if (allocated(this%boundaries)) deallocate(this%boundaries)
      allocate(this%boundaries(size(selectIndices)+1))
      ! set the boundaries based on the selection used
      this%boundaries(1) = periodValues(selectIndices(1))
      do iSel=1, size(selectIndices)
        this%boundaries(iSel+1) = periodValues(selectIndices(iSel)+1)
      end do
      ! set inital counter value to that of start year
      this%i = this%yearIds(simPerArg%yStart)
      this%nIds = size(selectIndices)
    case default
      ! this should be a string like "days since 1:2:3"
      if (.not. present(keepUnneededPeriods)) then
        call error_message('Cannot init aggregate period by name "', trim(units), &
                '", without flag keepUnneededPeriods.')
      end if
      ! check the units attribute and convert all periodValues to secondSince
      call check_time_unit(units, periodValues, this%refJulDate, inPeriod, this%secondsSince)
      this%isAveraged = .false.
      ! get correct selection of for the given simPerArg
      call get_period_indices(simPerArg, this%secondsSince, inPeriod, keepUnneededPeriods, &
              this%timeStep, this%i, selectIndices)
      this%nIds = size(selectIndices)
    end select

  end subroutine init_aggregateperiod

  subroutine increment_aggregateperiod(this, datetimeinfoArg)
    !< increment an aggregatePeriod counter based on the datetimeinfo argument
    !> aggregateperiod to increment
    class(aggregateperiod), intent(inout) :: this
    !> datetimeinfo to be used as reference
    class(datetimeinfo), intent(in) :: datetimeinfoArg

    if (this%isAveraged) then
      select case (this%timeStep)
      case(-1) ! day of year
        this%i = datetimeinfoArg%day
      case(-2) ! month of year
        this%i = datetimeinfoArg%month
      end select
    else
      select case (this%timeStep)
      case(-1) ! daily timestep
        if (datetimeinfoArg%is_new_day) then
          this%i = this%i + 1
        end if
      case(-2) ! monthly timestep
        if (datetimeinfoArg%is_new_month) then
          this%i = this%i + 1
        end if
      case(-3) ! yearly timestep
        if (datetimeinfoArg%is_new_year) then
          if (.not. this%isRegular) then
            ! this is hackish thing as it might be that this routine is called after the domainDateTime is incremented
            ! over the last valid this%years
            ! TODO: change this one the increment routines are called from the same spot in mhm_eval
            if (datetimeinfoArg%year <= ubound(this%yearIds, dim=1)) then
              this%i = this%yearIds(datetimeinfoArg%year)
            end if
          else
            this%i = this%i + 1
          end if
        end if
      end select
    end if

  end subroutine increment_aggregateperiod

  function get_unit_aggregateperiod(this) result(units)
    !< get unit attribute to be put to netcdf file
    use mo_julian, only: dec2date

    !> aggregateperiod containing the unit information
    class(aggregateperiod), intent(inout) :: this
    !> time unit descriptor
    character(256) :: units

    integer(i4) :: yy, mm, dd

    if (this%isAveraged) then
      select case (this%timeStep)
      case(-1) ! day of year
        units = 'day of year'
      case(-2) ! month of year
        units = 'month of year'
      end select
    else
      if (.not. this%isRegular) then
        write(units, '(A)') 'years'
      else
        ! this is the only CF conventions compatible attribute so far
        call dec2date(real(this%refJulDate, kind=dp), dd, mm, yy)
        write(units, '(A,I0.4,A,I0.2,A,I0.2)') 'days since ',yy, '-', mm, '-', dd
      end if
    end if

  end function get_unit_aggregateperiod

  function get_values_aggregateperiod(this) result(values)
    !< get unit attribute to be put to netcdf file

    use mo_constants, only: DaySecs

    !> aggregateperiod containing the value information
    class(aggregateperiod), intent(inout) :: this
    !> vector of period values (bounds for each period, thus size nIds+1)
    integer(i4), dimension(:), allocatable :: values

    integer(i4) :: iVal

    if (this%isAveraged) then
      values = [(iVal, iVal=1_i4, this%nIds + 1)]
    else
      if (.not. this%isRegular) then
        values = this%boundaries
      else
        values = int(this%secondsSince / DaySecs, kind=i4)
      end if
    end if

  end function get_values_aggregateperiod

  subroutine check_time_unit(units, periodValues, jRef, inPeriod, periodValuesSeconds)
    !< reads units string, converts periodValues to unit "seconds since date" (periodValuesSeconds, jRef) and also
    !< outputs the period covered by periodValues
    use mo_constants, only : DayHours, DaySecs
    use mo_julian, only : caldat, dec2date, julday
    use mo_kind, only : i8
    use mo_string_utils, only : divide_string

    !> time units descriptor
    character(*), intent(in) :: units
    !> periodValues for each Timestamp (center of interval/period)
    integer(i4), dimension(:), intent(in) :: periodValues
    !> reference julian date
    integer(i4), intent(out) :: jRef
    !> period covered by periodValues 
    class(period), intent(out), optional :: inPeriod
    !> the values in seconds and i8, thus seperate from periodValues
    integer(i8), dimension(:), allocatable, intent(out) :: periodValuesSeconds

    ! reference time
    integer(i4) :: yRef, dRef, mRef, hRef
    ! helper variables
    integer(i8) :: timeStepSeconds
    integer(i4) :: nTime
    integer(i4) :: hstart_int, hend_int
    ! dummies for netcdf attribute handling
    character(256), dimension(:), allocatable :: strArr, date, time

    ! units looks like "<unit> since YYYY-MM-DD[ HH:MM:SS]"
    ! split at space
    call divide_string(trim(units), ' ', strArr)

    ! determine reference time at '-' and convert to integer
    call divide_string(trim(strArr(3)), '-', date)
    read(date(1), *) yRef
    read(date(2), *) mRef
    read(date(3), *) dRef

    jRef = julday(dd = dRef, mm = mRef, yy = yRef)

    ! if existing also read in the time (only hour so far)
    hRef = 0
    if(size(strArr) > 3) then
      call divide_string(trim(strArr(4)), ':', time)
      read(time(1), *) hRef
    end if

    ! determine the step_size
    if (strArr(1) == 'days') then
      timeStepSeconds = int(DaySecs, kind=i8)
    else if (strArr(1) == 'hours') then
      timeStepSeconds = int(DaySecs / DayHours, kind=i8)
    else if (strArr(1) == 'minutes') then
      timeStepSeconds = int(DaySecs / DayHours / 60._dp, kind=i8)
    else if (strArr(1) == 'seconds') then
      timeStepSeconds = 1_i8
    else
      ! to suppress compiler warnings
      timeStepSeconds = -9999_i8
      call error_message('***ERROR: Please provide the input data in (days, hours, minutes, seconds) ', &
              'since YYYY-MM-DD[ HH:MM:SS] in the netcdf file. Found: ', trim(units))
    end if
    
    ! convert array from units since to seconds
    periodValuesSeconds = int(periodValues, kind=i8) * timeStepSeconds

    if (present(inPeriod)) then
      ! check for length of time vector, needs to be at least of length 2, otherwise step width check fails
      if (size(periodValues) < 2_i4) then
        call error_message('***ERROR: length of time dimension needs to be at least 2')
      end if
      ! compare the read period from ncfile to the period required
      ! convert julian second information back to date via conversion to float
      ! the 0.5_dp is for the different reference of fractional julian days, hours are truncated
      nTime = size(periodValuesSeconds)
      call dec2date(periodValuesSeconds(1) / DaySecs - 0.5_dp + jRef + hRef / 24._dp, inPeriod%dStart, inPeriod%mStart, &
              inPeriod%yStart, hstart_int)
      inPeriod%julStart = int(periodValuesSeconds(1) / DaySecs + jRef + hRef / 24._dp)
      call dec2date(periodValuesSeconds(nTime) / DaySecs - 0.5_dp + jRef + hRef / 24._dp, inPeriod%dEnd, inPeriod%mEnd, &
              inPeriod%yEnd, hend_int)
      inPeriod%julEnd = int(periodValuesSeconds(nTime) / DaySecs + jRef + hRef / 24._dp)
    end if

  end subroutine check_time_unit
    
  subroutine get_land_cover_period_indices(simPerArg, boundaries, keepUnneededPeriods, yearIds, selectIndices)
    !< compare periods defined by yearly boundary values to simulation period and retrieve vector of indices
    use mo_string_utils, only: compress
    use mo_message, only: error_message
    use mo_string_utils, only: num2str

    !> simulation period to be compared against
    class(period), intent(in) :: simPerArg
    !> year boundaries describing the aggregationperiod
    integer(i4), dimension(:), intent(in) :: boundaries
    !> flag describing whether to intersect boundaries with simPerArg
    logical, intent(in) :: keepUnneededPeriods
    !> vector of indices of the period for each simulation period year
    integer(i4), dimension(:), allocatable, intent(out), optional :: yearIds
    !> vector of indices of the periods to be selected
    integer(i4), dimension(:), allocatable, intent(out) :: selectIndices

    logical, dimension(size(boundaries) - 1) :: select_indices_mask
    logical, dimension(size(boundaries)) :: select_indices_temp

    integer(i4) :: select_index, iBoundary, LCyearStart, LCyearEnd

    select_index = 0_i4
    select_indices_mask = .false.
    if (present(yearIds)) then
      allocate(yearIds(simPerArg%yStart: simPerArg%yEnd))
    end if
    LCyearStart = simPerArg%yStart
    LCyearEnd = simPerArg%yEnd
    allocate(selectIndices(size(boundaries) - 1))

    ! set the correct indices to use
    do iBoundary=1, size(boundaries) - 1
      ! check for overlap ((StartA <= EndB) and (EndA >= StartB))
      ! https://stackoverflow.com/questions/325933/
      if ((boundaries(iBoundary) <= simPerArg%yend) .and. &
             ((boundaries(iBoundary+1)-1) >= simPerArg%ystart) .or. keepUnneededPeriods) then
        ! advance counter
        select_index = select_index + 1_i4
        ! select this iBoundary from dimension
        select_indices_mask(iBoundary) = .true.
        ! set the correct yearIds
        if (present(yearIds)) then
          yearIds(&
                  maxval([boundaries(iBoundary), LCyearStart]):&
                  minval([boundaries(iBoundary+1), LCyearEnd])) = select_index

        end if
      end if
    end do
    selectIndices = pack([(iBoundary, iBoundary=1, size(boundaries) - 1)], select_indices_mask)

    ! check if both start and end are covered
    ! this effectively selects the start boundaries (yStart)
    select_indices_temp = [select_indices_mask, .false.]
    if (minval(boundaries, mask=select_indices_temp) > simPerArg%ystart) then
      call error_message('The selected land cover periods ', &
              ' (', compress(trim(num2str(minval(boundaries, mask=select_indices_temp)))), &
              ') do not cover the beginning of the simulation period (', &
              compress(trim(num2str(simPerArg%ystart))), ').')
    end if
    ! this effectively selects the end boundaries (yEnd)
    select_indices_temp = [.false., select_indices_mask]
    if (maxval(boundaries, mask=select_indices_temp) < simPerArg%yend) then
      call error_message('The selected land cover periods ', &
              ' (', compress(trim(num2str(maxval(boundaries, mask=select_indices_temp)))), &
              ' ) do not cover the end of the simulation period (', &
              compress(trim(num2str(simPerArg%yend))), ').')
    end if

  end subroutine get_land_cover_period_indices

  subroutine get_period_indices(simPerArg, periodValues, inPeriod, keepUnneededPeriods, timeStep, i, selectIndices)
    !< \brief Extract time vector
    !< \details Extract time vector in unit julian hours and get supposed time step in hours

    use mo_constants, only : DaySecs, YearDays
    use mo_julian, only : caldat, julday
    use mo_kind, only : i8

    !> simulation period
    class(period), intent(in) :: simPerArg
    !> vector of date values [s]
    integer(i8), dimension(:), allocatable, intent(inout) :: periodValues
    !> period of input data which needs to be indexed/selected
    class(period), intent(in) :: inPeriod
    !> flag indicating whether to select only needed periods (defined by simulationperiod)
    logical, intent(in) :: keepUnneededPeriods
    !> time step of date values
    integer(i4), intent(out) :: timeStep
    !> start index for counter
    integer(i4), intent(out) :: i
    !> indices to select from inPeriod
    integer(i4), dimension(:), allocatable, intent(out) :: selectIndices


    real(dp), dimension(size(periodValues)-1) :: periodValuesDiff
    integer(i4) :: ncJulSta1, dd, nTime, iInd, nSel, iErr
    integer(i4) :: mmcalstart, mmcalend, yycalstart, yycalend
    integer(i4) :: mmncstart, yyncstart

    ! default in case of use keepUnneededPeriods: use all indices of inPeriod
    nTime = size(periodValues)
    selectIndices = [(iInd, iInd=1, nTime-1)]
    periodValuesDiff = real(abs(periodValues(2 : nTime) - periodValues(1 : nTime - 1)), kind=dp)

    ! prepare the selection and check for required time_step
    ! sanity check: difference must be 1 day
    if (all((periodValuesDiff / DaySecs - 1._dp) <= 1.e-6)) then
      ncJulSta1 = inPeriod%julStart
      ! in case of use keepUnneededPeriods: this is the first valid index for simPerArg
      i = simPerArg%julStart - ncJulSta1 + 1_i4
      nSel = simPerArg%julEnd - simPerArg%julStart + 1_i4
      timeStep = -1
    ! sanity check: difference must be between 28 and 31 days
    elseif (all((periodValuesDiff / DaySecs) <= 31._dp) .and. &
            all((periodValuesDiff / DaySecs) >= 28._dp)) then
      call caldat(simPerArg%julStart, dd, mmcalstart, yycalstart)
      call caldat(inPeriod%julStart, dd, mmncstart, yyncstart)
      ! monthly timesteps are usually set by month end, so for beginning, we need 1st of month
      ncJulSta1 = julday(1, mmncstart, yyncstart)
      call caldat(simPerArg%julEnd, dd, mmcalend, yycalend)
      ! in case of use keepUnneededPeriods: this is the first valid index for simPerArg
      i = (yycalstart * 12 + mmcalstart) - (yyncstart * 12 + mmncstart) + 1_i4
      nSel = (yycalend * 12 + mmcalend) - (yycalstart * 12 + mmcalstart) + 1_i4
      timeStep = -2  ! monthly time step
    ! difference must be between 365 and 366 days
    elseif (all((periodValuesDiff / DaySecs) <= (YearDays + 1._dp)) .and. &
            all((periodValuesDiff / DaySecs) >= YearDays)) then
      call caldat(simPerArg%julStart, dd, mmcalstart, yycalstart)
      call caldat(inPeriod%julStart, dd, mmncstart, yyncstart)
      ! yearly timesteps are usually set by year end, so for beginning, we need 1st of year
      ncJulSta1 = julday(1, 1, yyncstart)
      call caldat(simPerArg%julEnd, dd, mmcalend, yycalend)
      ! in case of use keepUnneededPeriods: this is the first valid index for simPerArg
      i = yycalstart - yyncstart + 1_i4
      nSel = yycalend - yycalstart + 1_i4
      timeStep = -3  ! yearly time step
    ! difference must be 1 hour
    elseif (all((periodValuesDiff / 3600._dp - 1._dp) <= 1.e-6)) then
      ncJulSta1 = inPeriod%julStart
      i = (simPerArg%julStart - ncJulSta1) * 24_i4 + 1_i4 ! convert to hours; always starts at one
      nSel = (simPerArg%julEnd - simPerArg%julStart + 1_i4) * 24_i4 ! convert to hours
      timeStep = -4  ! hourly time step
    else
      call message('First time step differences in file [days]:')
      do iErr=1, max(3_i4, size(periodValuesDiff))
        call message(trim(num2str(iErr)), ' :', trim(num2str(periodValuesDiff(iErr) / DaySecs)))
      end do
      call error_message('***ERROR: time step cannot be inferred because steps are not equal over all times ', &
              'and/or do not conform to any common timestep.' )
    end if
    if (.not. keepUnneededPeriods) then
      ! we select only the relevant slice of dates
      selectIndices = selectIndices(i: i+nSel-1)
      ! apply it to values stored
      periodValues = periodValues(i: i+nSel)
      ! ... and thus start with index 1
      i = 1_i4
    end if

    ! Check if time steps in file cover simulation period
    if (.not. ((ncJulSta1 <= simPerArg%julStart) .AND. (inPeriod%julEnd >= simPerArg%julEnd))) then
      call error_message('***ERROR: time period of input data (', &
              trim(num2str(ncJulSta1)), trim(num2str(inPeriod%julEnd)), ') is not matching modelling period (',&
              trim(num2str(simPerArg%julStart)), trim(num2str(simPerArg%julEnd)), ').')
    end if

  end subroutine get_period_indices

END MODULE mo_common_datetime_type
