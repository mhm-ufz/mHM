!> \file nml_parameter_lakes_1.f90
!> \copydoc nml_lakes_1

!> \brief Lakes - Case 1
!> \details Parameters for the mHM lake module (mLM) controlling spill over the dam crest and the
!! seepage and baseflow from the lake groundwater store. Each parameter is given per lake
!! per domain as the four components value, optimize-flag, min, max (min/max are the
!! calibration bounds). Reservoir release-regulation parameters are not included at this
!! stage (minimal natural-lake scope).
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_lakes_1
  use nml_helper, only: &
    nml_file_t, &
    nml_line_buffer, &
    NML_OK, &
    NML_ERR_FILE_NOT_FOUND, &
    NML_ERR_OPEN, &
    NML_ERR_NOT_OPEN, &
    NML_ERR_NML_NOT_FOUND, &
    NML_ERR_READ, &
    NML_ERR_CLOSE, &
    NML_ERR_REQUIRED, &
    NML_ERR_ENUM, &
    NML_ERR_BOUNDS, &
    NML_ERR_NOT_SET, &
    NML_ERR_INVALID_NAME, &
    NML_ERR_INVALID_INDEX, &
    idx_check, &
    to_lower, &
    n_domains__default, &
    n_lakes__default, &
    NML_ERR_PARTLY_SET
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp
  use mo_parameter_types, only: parameter_t

  implicit none

  !> \class nml_lakes_1_t
  !> \brief Lakes - Case 1
  !> \details Parameters for the mHM lake module (mLM) controlling spill over the dam crest and the
  !! seepage and baseflow from the lake groundwater store. Each parameter is given per lake
  !! per domain as the four components value, optimize-flag, min, max (min/max are the
  !! calibration bounds). Reservoir release-regulation parameters are not included at this
  !! stage (minimal natural-lake scope).
  type, public :: nml_lakes_1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    integer :: n_lakes = n_lakes__default !< runtime dimension for n_lakes
    type(parameter_t), allocatable, dimension(:, :) :: spill_coefficient_of_discharge !< Spill coefficient of discharge
    type(parameter_t), allocatable, dimension(:, :) :: spill_exponent !< Spill exponent
    type(parameter_t), allocatable, dimension(:, :) :: longterm_baseflow !< Long-term baseflow
    type(parameter_t), allocatable, dimension(:, :) :: exponent_for_percolation !< Exponent for percolation
  contains
    procedure :: init => nml_lakes_1_init
    procedure :: init_type => nml_lakes_1_init_type
    procedure :: set_dims => nml_lakes_1_set_dims
    procedure :: from_file => nml_lakes_1_from_file
    procedure :: set => nml_lakes_1_set
    procedure :: is_set => nml_lakes_1_is_set
    procedure :: is_valid => nml_lakes_1_is_valid
  end type nml_lakes_1_t

contains

  !> \brief Initialize defaults and sentinels for lakes_1
  integer function nml_lakes_1_init(this, errmsg) result(status)
    class(nml_lakes_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      spill_coefficient_of_discharge=this%spill_coefficient_of_discharge, &
      spill_exponent=this%spill_exponent, &
      longterm_baseflow=this%longterm_baseflow, &
      exponent_for_percolation=this%exponent_for_percolation, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_lakes_1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_lakes_1_init_type(this, &
    spill_coefficient_of_discharge, &
    spill_exponent, &
    longterm_baseflow, &
    exponent_for_percolation, &
    errmsg) result(status)
    class(nml_lakes_1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), dimension(:, :), allocatable, intent(inout), optional :: spill_coefficient_of_discharge !< Spill coefficient of discharge
    type(parameter_t), dimension(:, :), allocatable, intent(inout), optional :: spill_exponent !< Spill exponent
    type(parameter_t), dimension(:, :), allocatable, intent(inout), optional :: longterm_baseflow !< Long-term baseflow
    type(parameter_t), dimension(:, :), allocatable, intent(inout), optional :: exponent_for_percolation !< Exponent for percolation
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(spill_coefficient_of_discharge)) then
      if (allocated(spill_coefficient_of_discharge)) deallocate(spill_coefficient_of_discharge)
      allocate(spill_coefficient_of_discharge(this%n_domains, this%n_lakes))
      spill_coefficient_of_discharge%value = ieee_value(spill_coefficient_of_discharge%value, ieee_quiet_nan) ! sentinel for derived component value
      spill_coefficient_of_discharge%optimize = .false.
      spill_coefficient_of_discharge%min = ieee_value(spill_coefficient_of_discharge%min, ieee_quiet_nan) ! sentinel for derived component min
      spill_coefficient_of_discharge%max = ieee_value(spill_coefficient_of_discharge%max, ieee_quiet_nan) ! sentinel for derived component max
      spill_coefficient_of_discharge%min = 1.4_dp
      spill_coefficient_of_discharge%max = 3.0_dp
    end if
    if (present(spill_exponent)) then
      if (allocated(spill_exponent)) deallocate(spill_exponent)
      allocate(spill_exponent(this%n_domains, this%n_lakes))
      spill_exponent%value = ieee_value(spill_exponent%value, ieee_quiet_nan) ! sentinel for derived component value
      spill_exponent%optimize = .false.
      spill_exponent%min = ieee_value(spill_exponent%min, ieee_quiet_nan) ! sentinel for derived component min
      spill_exponent%max = ieee_value(spill_exponent%max, ieee_quiet_nan) ! sentinel for derived component max
      spill_exponent%min = 0.5_dp
      spill_exponent%max = 2.5_dp
    end if
    if (present(longterm_baseflow)) then
      if (allocated(longterm_baseflow)) deallocate(longterm_baseflow)
      allocate(longterm_baseflow(this%n_domains, this%n_lakes))
      longterm_baseflow%value = ieee_value(longterm_baseflow%value, ieee_quiet_nan) ! sentinel for derived component value
      longterm_baseflow%optimize = .false.
      longterm_baseflow%min = ieee_value(longterm_baseflow%min, ieee_quiet_nan) ! sentinel for derived component min
      longterm_baseflow%max = ieee_value(longterm_baseflow%max, ieee_quiet_nan) ! sentinel for derived component max
      longterm_baseflow%min = 0.0001_dp
      longterm_baseflow%max = 5.0_dp
    end if
    if (present(exponent_for_percolation)) then
      if (allocated(exponent_for_percolation)) deallocate(exponent_for_percolation)
      allocate(exponent_for_percolation(this%n_domains, this%n_lakes))
      exponent_for_percolation%value = ieee_value(exponent_for_percolation%value, ieee_quiet_nan) ! sentinel for derived component value
      exponent_for_percolation%optimize = .false.
      exponent_for_percolation%min = ieee_value(exponent_for_percolation%min, ieee_quiet_nan) ! sentinel for derived component min
      exponent_for_percolation%max = ieee_value(exponent_for_percolation%max, ieee_quiet_nan) ! sentinel for derived component max
      exponent_for_percolation%min = 8.0_dp
      exponent_for_percolation%max = 15.0_dp
    end if
  end function nml_lakes_1_init_type

  !> \brief Reset runtime dimensions for lakes_1
  integer function nml_lakes_1_set_dims(this, &
    n_domains, &
    n_lakes, &
    errmsg) result(status)
    class(nml_lakes_1_t), intent(inout) :: this !< namelist instance
    integer, intent(in), optional :: n_domains !< runtime dimension override for n_domains
    integer, intent(in), optional :: n_lakes !< runtime dimension override for n_lakes
    integer :: candidate__n_domains
    integer :: candidate__n_lakes
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(n_domains)) then
      candidate__n_domains = n_domains
    else
      candidate__n_domains = n_domains__default
    end if
    if (candidate__n_domains <= 0) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 'n_domains' must be positive"
      return
    end if
    if (present(n_lakes)) then
      candidate__n_lakes = n_lakes
    else
      candidate__n_lakes = n_lakes__default
    end if
    if (candidate__n_lakes <= 0) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 'n_lakes' must be positive"
      return
    end if
    this%n_domains = candidate__n_domains
    this%n_lakes = candidate__n_lakes

    ! deallocate runtime-sized fields; init/set/from_file allocate them again
    if (allocated(this%spill_coefficient_of_discharge)) deallocate(this%spill_coefficient_of_discharge)
    if (allocated(this%spill_exponent)) deallocate(this%spill_exponent)
    if (allocated(this%longterm_baseflow)) deallocate(this%longterm_baseflow)
    if (allocated(this%exponent_for_percolation)) deallocate(this%exponent_for_percolation)
    this%is_configured = .false.
  end function nml_lakes_1_set_dims


  !> \brief Read lakes_1 namelist from file
  integer function nml_lakes_1_from_file(this, file, errmsg) result(status)
    class(nml_lakes_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t), allocatable, dimension(:, :) :: spill_coefficient_of_discharge
    type(parameter_t), allocatable, dimension(:, :) :: spill_exponent
    type(parameter_t), allocatable, dimension(:, :) :: longterm_baseflow
    type(parameter_t), allocatable, dimension(:, :) :: exponent_for_percolation
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /lakes_1/ &
      spill_coefficient_of_discharge, &
      spill_exponent, &
      longterm_baseflow, &
      exponent_for_percolation

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(spill_coefficient_of_discharge)) deallocate(spill_coefficient_of_discharge)
    allocate(spill_coefficient_of_discharge(this%n_domains, this%n_lakes))
    if (allocated(spill_exponent)) deallocate(spill_exponent)
    allocate(spill_exponent(this%n_domains, this%n_lakes))
    if (allocated(longterm_baseflow)) deallocate(longterm_baseflow)
    allocate(longterm_baseflow(this%n_domains, this%n_lakes))
    if (allocated(exponent_for_percolation)) deallocate(exponent_for_percolation)
    allocate(exponent_for_percolation(this%n_domains, this%n_lakes))
    spill_coefficient_of_discharge = this%spill_coefficient_of_discharge
    spill_exponent = this%spill_exponent
    longterm_baseflow = this%longterm_baseflow
    exponent_for_percolation = this%exponent_for_percolation

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("lakes_1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=lakes_1, iostat=iostat, iomsg=iomsg)
    if (iostat /= 0) then
      status = NML_ERR_READ
      if (present(errmsg)) errmsg = trim(iomsg)
      close_status = nml%close()
      return
    end if
    close_status = nml%close(errmsg=errmsg)
    if (close_status /= NML_OK) then
      status = close_status
      return
    end if

    ! assign values
    this%spill_coefficient_of_discharge = spill_coefficient_of_discharge
    this%spill_exponent = spill_exponent
    this%longterm_baseflow = longterm_baseflow
    this%exponent_for_percolation = exponent_for_percolation

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_lakes_1_from_file

  !> \brief Set lakes_1 values
  integer function nml_lakes_1_set(this, &
    spill_coefficient_of_discharge, &
    spill_exponent, &
    longterm_baseflow, &
    exponent_for_percolation, &
    errmsg) result(status)

    class(nml_lakes_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), dimension(:, :), intent(in) :: spill_coefficient_of_discharge !< Spill coefficient of discharge
    type(parameter_t), dimension(:, :), intent(in) :: spill_exponent !< Spill exponent
    type(parameter_t), dimension(:, :), intent(in) :: longterm_baseflow !< Long-term baseflow
    type(parameter_t), dimension(:, :), intent(in) :: exponent_for_percolation !< Exponent for percolation
    integer :: &
      lb__1, &
      lb__2, &
      ub__1, &
      ub__2

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    if (size(spill_coefficient_of_discharge, 1) > size(this%spill_coefficient_of_discharge, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'spill_coefficient_of_discharge'"
      return
    end if
    lb__1 = lbound(this%spill_coefficient_of_discharge, 1)
    ub__1 = lb__1 + size(spill_coefficient_of_discharge, 1) - 1
    if (size(spill_coefficient_of_discharge, 2) > size(this%spill_coefficient_of_discharge, 2)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'spill_coefficient_of_discharge'"
      return
    end if
    lb__2 = lbound(this%spill_coefficient_of_discharge, 2)
    ub__2 = lb__2 + size(spill_coefficient_of_discharge, 2) - 1
    this%spill_coefficient_of_discharge(lb__1:ub__1, lb__2:ub__2) = spill_coefficient_of_discharge
    if (size(spill_exponent, 1) > size(this%spill_exponent, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'spill_exponent'"
      return
    end if
    lb__1 = lbound(this%spill_exponent, 1)
    ub__1 = lb__1 + size(spill_exponent, 1) - 1
    if (size(spill_exponent, 2) > size(this%spill_exponent, 2)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'spill_exponent'"
      return
    end if
    lb__2 = lbound(this%spill_exponent, 2)
    ub__2 = lb__2 + size(spill_exponent, 2) - 1
    this%spill_exponent(lb__1:ub__1, lb__2:ub__2) = spill_exponent
    if (size(longterm_baseflow, 1) > size(this%longterm_baseflow, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'longterm_baseflow'"
      return
    end if
    lb__1 = lbound(this%longterm_baseflow, 1)
    ub__1 = lb__1 + size(longterm_baseflow, 1) - 1
    if (size(longterm_baseflow, 2) > size(this%longterm_baseflow, 2)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'longterm_baseflow'"
      return
    end if
    lb__2 = lbound(this%longterm_baseflow, 2)
    ub__2 = lb__2 + size(longterm_baseflow, 2) - 1
    this%longterm_baseflow(lb__1:ub__1, lb__2:ub__2) = longterm_baseflow
    if (size(exponent_for_percolation, 1) > size(this%exponent_for_percolation, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'exponent_for_percolation'"
      return
    end if
    lb__1 = lbound(this%exponent_for_percolation, 1)
    ub__1 = lb__1 + size(exponent_for_percolation, 1) - 1
    if (size(exponent_for_percolation, 2) > size(this%exponent_for_percolation, 2)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'exponent_for_percolation'"
      return
    end if
    lb__2 = lbound(this%exponent_for_percolation, 2)
    ub__2 = lb__2 + size(exponent_for_percolation, 2) - 1
    this%exponent_for_percolation(lb__1:ub__1, lb__2:ub__2) = exponent_for_percolation

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_lakes_1_set

  !> \brief Check whether a namelist value was set
  integer function nml_lakes_1_is_set(this, name, idx, errmsg) result(status)
    class(nml_lakes_1_t), intent(in) :: this !< namelist instance
    character(len=*), intent(in) :: name !< field name
    integer, intent(in), optional :: idx(:) !< optional field index values
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if
    select case (to_lower(trim(name)))
    case ("spill_coefficient_of_discharge%value")
      if (.not. allocated(this%spill_coefficient_of_discharge)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_coefficient_of_discharge), ubound(this%spill_coefficient_of_discharge), &
          "spill_coefficient_of_discharge", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%spill_coefficient_of_discharge(idx(1), idx(2))%value)) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%spill_coefficient_of_discharge%value))) status = NML_ERR_NOT_SET
      end if
    case ("spill_coefficient_of_discharge%optimize")
      if (.not. allocated(this%spill_coefficient_of_discharge)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_coefficient_of_discharge), ubound(this%spill_coefficient_of_discharge), &
          "spill_coefficient_of_discharge", errmsg)
        if (status /= NML_OK) return
      end if
    case ("spill_coefficient_of_discharge%min")
      if (.not. allocated(this%spill_coefficient_of_discharge)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_coefficient_of_discharge), ubound(this%spill_coefficient_of_discharge), &
          "spill_coefficient_of_discharge", errmsg)
        if (status /= NML_OK) return
      end if
    case ("spill_coefficient_of_discharge%max")
      if (.not. allocated(this%spill_coefficient_of_discharge)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_coefficient_of_discharge), ubound(this%spill_coefficient_of_discharge), &
          "spill_coefficient_of_discharge", errmsg)
        if (status /= NML_OK) return
      end if
    case ("spill_coefficient_of_discharge")
      if (.not. allocated(this%spill_coefficient_of_discharge)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_coefficient_of_discharge), ubound(this%spill_coefficient_of_discharge), &
          "spill_coefficient_of_discharge", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%spill_coefficient_of_discharge(idx(1), idx(2))%value)) then
          status = NML_ERR_NOT_SET
        end if
      else
        if (all(ieee_is_nan(this%spill_coefficient_of_discharge%value))) then
          status = NML_ERR_NOT_SET
        else if (any(ieee_is_nan(this%spill_coefficient_of_discharge%value))) then
          status = NML_ERR_PARTLY_SET
        end if
      end if
    case ("spill_exponent%value")
      if (.not. allocated(this%spill_exponent)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_exponent), ubound(this%spill_exponent), &
          "spill_exponent", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%spill_exponent(idx(1), idx(2))%value)) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%spill_exponent%value))) status = NML_ERR_NOT_SET
      end if
    case ("spill_exponent%optimize")
      if (.not. allocated(this%spill_exponent)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_exponent), ubound(this%spill_exponent), &
          "spill_exponent", errmsg)
        if (status /= NML_OK) return
      end if
    case ("spill_exponent%min")
      if (.not. allocated(this%spill_exponent)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_exponent), ubound(this%spill_exponent), &
          "spill_exponent", errmsg)
        if (status /= NML_OK) return
      end if
    case ("spill_exponent%max")
      if (.not. allocated(this%spill_exponent)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_exponent), ubound(this%spill_exponent), &
          "spill_exponent", errmsg)
        if (status /= NML_OK) return
      end if
    case ("spill_exponent")
      if (.not. allocated(this%spill_exponent)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%spill_exponent), ubound(this%spill_exponent), &
          "spill_exponent", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%spill_exponent(idx(1), idx(2))%value)) then
          status = NML_ERR_NOT_SET
        end if
      else
        if (all(ieee_is_nan(this%spill_exponent%value))) then
          status = NML_ERR_NOT_SET
        else if (any(ieee_is_nan(this%spill_exponent%value))) then
          status = NML_ERR_PARTLY_SET
        end if
      end if
    case ("longterm_baseflow%value")
      if (.not. allocated(this%longterm_baseflow)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%longterm_baseflow), ubound(this%longterm_baseflow), &
          "longterm_baseflow", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%longterm_baseflow(idx(1), idx(2))%value)) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%longterm_baseflow%value))) status = NML_ERR_NOT_SET
      end if
    case ("longterm_baseflow%optimize")
      if (.not. allocated(this%longterm_baseflow)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%longterm_baseflow), ubound(this%longterm_baseflow), &
          "longterm_baseflow", errmsg)
        if (status /= NML_OK) return
      end if
    case ("longterm_baseflow%min")
      if (.not. allocated(this%longterm_baseflow)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%longterm_baseflow), ubound(this%longterm_baseflow), &
          "longterm_baseflow", errmsg)
        if (status /= NML_OK) return
      end if
    case ("longterm_baseflow%max")
      if (.not. allocated(this%longterm_baseflow)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%longterm_baseflow), ubound(this%longterm_baseflow), &
          "longterm_baseflow", errmsg)
        if (status /= NML_OK) return
      end if
    case ("longterm_baseflow")
      if (.not. allocated(this%longterm_baseflow)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%longterm_baseflow), ubound(this%longterm_baseflow), &
          "longterm_baseflow", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%longterm_baseflow(idx(1), idx(2))%value)) then
          status = NML_ERR_NOT_SET
        end if
      else
        if (all(ieee_is_nan(this%longterm_baseflow%value))) then
          status = NML_ERR_NOT_SET
        else if (any(ieee_is_nan(this%longterm_baseflow%value))) then
          status = NML_ERR_PARTLY_SET
        end if
      end if
    case ("exponent_for_percolation%value")
      if (.not. allocated(this%exponent_for_percolation)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%exponent_for_percolation), ubound(this%exponent_for_percolation), &
          "exponent_for_percolation", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%exponent_for_percolation(idx(1), idx(2))%value)) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%exponent_for_percolation%value))) status = NML_ERR_NOT_SET
      end if
    case ("exponent_for_percolation%optimize")
      if (.not. allocated(this%exponent_for_percolation)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%exponent_for_percolation), ubound(this%exponent_for_percolation), &
          "exponent_for_percolation", errmsg)
        if (status /= NML_OK) return
      end if
    case ("exponent_for_percolation%min")
      if (.not. allocated(this%exponent_for_percolation)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%exponent_for_percolation), ubound(this%exponent_for_percolation), &
          "exponent_for_percolation", errmsg)
        if (status /= NML_OK) return
      end if
    case ("exponent_for_percolation%max")
      if (.not. allocated(this%exponent_for_percolation)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%exponent_for_percolation), ubound(this%exponent_for_percolation), &
          "exponent_for_percolation", errmsg)
        if (status /= NML_OK) return
      end if
    case ("exponent_for_percolation")
      if (.not. allocated(this%exponent_for_percolation)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%exponent_for_percolation), ubound(this%exponent_for_percolation), &
          "exponent_for_percolation", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%exponent_for_percolation(idx(1), idx(2))%value)) then
          status = NML_ERR_NOT_SET
        end if
      else
        if (all(ieee_is_nan(this%exponent_for_percolation%value))) then
          status = NML_ERR_NOT_SET
        else if (any(ieee_is_nan(this%exponent_for_percolation%value))) then
          status = NML_ERR_PARTLY_SET
        end if
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_lakes_1_is_set

  !> \brief Validate required values and constraints
  integer function nml_lakes_1_is_valid(this, errmsg) result(status)
    class(nml_lakes_1_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! required parameters
    istat = this%is_set("spill_coefficient_of_discharge", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: spill_coefficient_of_discharge"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("spill_exponent", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: spill_exponent"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("longterm_baseflow", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: longterm_baseflow"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("exponent_for_percolation", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: exponent_for_percolation"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_lakes_1_is_valid

end module nml_lakes_1
