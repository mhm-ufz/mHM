!> \file nml_parameter_lakes_1.f90
!> \copydoc nml_lakes_1

!> \brief Lakes - Case 1
!> \details Parameters for the mHM lake module (mLM) controlling spill over the dam crest and the
!! seepage and baseflow from the lake groundwater store. All parameters are global (pan-domain, shared across all lakes, like the routing
!! parameters), so each is a single scalar given as the four components value, optimize-flag,
!! min, max (min/max are the calibration bounds).
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
    n_geo_units__default, &
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
  !! seepage and baseflow from the lake groundwater store. All parameters are global (pan-domain, shared across all lakes, like the routing
  !! parameters), so each is a single scalar given as the four components value, optimize-flag,
  !! min, max (min/max are the calibration bounds).
  type, public :: nml_lakes_1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_geo_units = n_geo_units__default !< runtime dimension for n_geo_units
    type(parameter_t) :: spill_coefficient_of_discharge !< Spill coefficient of discharge
    type(parameter_t) :: spill_exponent !< Spill exponent
    type(parameter_t) :: b !< Effective lakebed thickness
    type(parameter_t) :: ptf_ks_constant !< Lakebed PTF constant for saturated hydraulic conductivity
    type(parameter_t) :: ptf_ks_sand !< Lakebed PTF sand multiplier for saturated hydraulic conductivity
    type(parameter_t) :: ptf_ks_clay !< Lakebed PTF clay multiplier for saturated hydraulic conductivity
    type(parameter_t), allocatable, dimension(:) :: baseflow_recession !< Lake baseflow recession time for each geological unit [d]
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
      b=this%b, &
      ptf_ks_constant=this%ptf_ks_constant, &
      ptf_ks_sand=this%ptf_ks_sand, &
      ptf_ks_clay=this%ptf_ks_clay, &
      baseflow_recession=this%baseflow_recession, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_lakes_1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_lakes_1_init_type(this, &
    spill_coefficient_of_discharge, &
    spill_exponent, &
    b, &
    ptf_ks_constant, &
    ptf_ks_sand, &
    ptf_ks_clay, &
    baseflow_recession, &
    errmsg) result(status)
    class(nml_lakes_1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: spill_coefficient_of_discharge !< Spill coefficient of discharge
    type(parameter_t), intent(inout), optional :: spill_exponent !< Spill exponent
    type(parameter_t), intent(inout), optional :: b !< Effective lakebed thickness
    type(parameter_t), intent(inout), optional :: ptf_ks_constant !< Lakebed PTF constant for saturated hydraulic conductivity
    type(parameter_t), intent(inout), optional :: ptf_ks_sand !< Lakebed PTF sand multiplier for saturated hydraulic conductivity
    type(parameter_t), intent(inout), optional :: ptf_ks_clay !< Lakebed PTF clay multiplier for saturated hydraulic conductivity
    type(parameter_t), dimension(:), allocatable, intent(inout), optional :: baseflow_recession !< Lake baseflow recession time for each geological unit [d]
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(spill_coefficient_of_discharge)) then
      spill_coefficient_of_discharge%value = ieee_value(spill_coefficient_of_discharge%value, ieee_quiet_nan) ! sentinel for derived component value
      spill_coefficient_of_discharge%optimize = .false.
      spill_coefficient_of_discharge%min = ieee_value(spill_coefficient_of_discharge%min, ieee_quiet_nan) ! sentinel for derived component min
      spill_coefficient_of_discharge%max = ieee_value(spill_coefficient_of_discharge%max, ieee_quiet_nan) ! sentinel for derived component max
      spill_coefficient_of_discharge%min = 1.4_dp
      spill_coefficient_of_discharge%max = 3.0_dp
    end if
    if (present(spill_exponent)) then
      spill_exponent%value = ieee_value(spill_exponent%value, ieee_quiet_nan) ! sentinel for derived component value
      spill_exponent%optimize = .false.
      spill_exponent%min = ieee_value(spill_exponent%min, ieee_quiet_nan) ! sentinel for derived component min
      spill_exponent%max = ieee_value(spill_exponent%max, ieee_quiet_nan) ! sentinel for derived component max
      spill_exponent%min = 0.5_dp
      spill_exponent%max = 2.5_dp
    end if
    if (present(b)) then
      b%value = ieee_value(b%value, ieee_quiet_nan) ! sentinel for derived component value
      b%optimize = .false.
      b%min = ieee_value(b%min, ieee_quiet_nan) ! sentinel for derived component min
      b%max = ieee_value(b%max, ieee_quiet_nan) ! sentinel for derived component max
      b%min = 0.1_dp
      b%max = 10.0_dp
    end if
    if (present(ptf_ks_constant)) then
      ptf_ks_constant%value = ieee_value(ptf_ks_constant%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_ks_constant%optimize = .false.
      ptf_ks_constant%min = ieee_value(ptf_ks_constant%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_ks_constant%max = ieee_value(ptf_ks_constant%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_ks_constant%min = -1.2_dp
      ptf_ks_constant%max = -0.285_dp
    end if
    if (present(ptf_ks_sand)) then
      ptf_ks_sand%value = ieee_value(ptf_ks_sand%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_ks_sand%optimize = .false.
      ptf_ks_sand%min = ieee_value(ptf_ks_sand%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_ks_sand%max = ieee_value(ptf_ks_sand%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_ks_sand%min = 0.006_dp
      ptf_ks_sand%max = 0.026_dp
    end if
    if (present(ptf_ks_clay)) then
      ptf_ks_clay%value = ieee_value(ptf_ks_clay%value, ieee_quiet_nan) ! sentinel for derived component value
      ptf_ks_clay%optimize = .false.
      ptf_ks_clay%min = ieee_value(ptf_ks_clay%min, ieee_quiet_nan) ! sentinel for derived component min
      ptf_ks_clay%max = ieee_value(ptf_ks_clay%max, ieee_quiet_nan) ! sentinel for derived component max
      ptf_ks_clay%min = 0.003_dp
      ptf_ks_clay%max = 0.013_dp
    end if
    if (present(baseflow_recession)) then
      if (allocated(baseflow_recession)) deallocate(baseflow_recession)
      allocate(baseflow_recession(this%n_geo_units))
      baseflow_recession%value = ieee_value(baseflow_recession%value, ieee_quiet_nan) ! sentinel for derived component value
      baseflow_recession%optimize = .false.
      baseflow_recession%min = ieee_value(baseflow_recession%min, ieee_quiet_nan) ! sentinel for derived component min
      baseflow_recession%max = ieee_value(baseflow_recession%max, ieee_quiet_nan) ! sentinel for derived component max
      baseflow_recession%min = 1.0_dp
      baseflow_recession%max = 1000.0_dp
    end if
  end function nml_lakes_1_init_type

  !> \brief Reset runtime dimensions for lakes_1
  integer function nml_lakes_1_set_dims(this, &
    n_geo_units, &
    errmsg) result(status)
    class(nml_lakes_1_t), intent(inout) :: this !< namelist instance
    integer, intent(in), optional :: n_geo_units !< runtime dimension override for n_geo_units
    integer :: candidate__n_geo_units
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(n_geo_units)) then
      candidate__n_geo_units = n_geo_units
    else
      candidate__n_geo_units = n_geo_units__default
    end if
    if (candidate__n_geo_units <= 0) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 'n_geo_units' must be positive"
      return
    end if
    this%n_geo_units = candidate__n_geo_units

    ! deallocate runtime-sized fields; init/set/from_file allocate them again
    if (allocated(this%baseflow_recession)) deallocate(this%baseflow_recession)
    this%is_configured = .false.
  end function nml_lakes_1_set_dims


  !> \brief Read lakes_1 namelist from file
  integer function nml_lakes_1_from_file(this, file, errmsg) result(status)
    class(nml_lakes_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: spill_coefficient_of_discharge
    type(parameter_t) :: spill_exponent
    type(parameter_t) :: b
    type(parameter_t) :: ptf_ks_constant
    type(parameter_t) :: ptf_ks_sand
    type(parameter_t) :: ptf_ks_clay
    type(parameter_t), allocatable, dimension(:) :: baseflow_recession
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /lakes_1/ &
      spill_coefficient_of_discharge, &
      spill_exponent, &
      b, &
      ptf_ks_constant, &
      ptf_ks_sand, &
      ptf_ks_clay, &
      baseflow_recession

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(baseflow_recession)) deallocate(baseflow_recession)
    allocate(baseflow_recession(this%n_geo_units))
    spill_coefficient_of_discharge = this%spill_coefficient_of_discharge
    spill_exponent = this%spill_exponent
    b = this%b
    ptf_ks_constant = this%ptf_ks_constant
    ptf_ks_sand = this%ptf_ks_sand
    ptf_ks_clay = this%ptf_ks_clay
    baseflow_recession = this%baseflow_recession

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
    this%b = b
    this%ptf_ks_constant = ptf_ks_constant
    this%ptf_ks_sand = ptf_ks_sand
    this%ptf_ks_clay = ptf_ks_clay
    this%baseflow_recession = baseflow_recession

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_lakes_1_from_file

  !> \brief Set lakes_1 values
  integer function nml_lakes_1_set(this, &
    spill_coefficient_of_discharge, &
    spill_exponent, &
    b, &
    ptf_ks_constant, &
    ptf_ks_sand, &
    ptf_ks_clay, &
    baseflow_recession, &
    errmsg) result(status)

    class(nml_lakes_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: spill_coefficient_of_discharge !< Spill coefficient of discharge
    type(parameter_t), intent(in) :: spill_exponent !< Spill exponent
    type(parameter_t), intent(in) :: b !< Effective lakebed thickness
    type(parameter_t), intent(in) :: ptf_ks_constant !< Lakebed PTF constant for saturated hydraulic conductivity
    type(parameter_t), intent(in) :: ptf_ks_sand !< Lakebed PTF sand multiplier for saturated hydraulic conductivity
    type(parameter_t), intent(in) :: ptf_ks_clay !< Lakebed PTF clay multiplier for saturated hydraulic conductivity
    type(parameter_t), dimension(:), intent(in) :: baseflow_recession !< Lake baseflow recession time for each geological unit [d]
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%spill_coefficient_of_discharge = spill_coefficient_of_discharge
    this%spill_exponent = spill_exponent
    this%b = b
    this%ptf_ks_constant = ptf_ks_constant
    this%ptf_ks_sand = ptf_ks_sand
    this%ptf_ks_clay = ptf_ks_clay
    if (size(baseflow_recession, 1) > size(this%baseflow_recession, 1)) then
      status = NML_ERR_INVALID_INDEX
      if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'baseflow_recession'"
      return
    end if
    lb__1 = lbound(this%baseflow_recession, 1)
    ub__1 = lb__1 + size(baseflow_recession, 1) - 1
    this%baseflow_recession(lb__1:ub__1) = baseflow_recession

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
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_coefficient_of_discharge'"
        return
      end if
      if (ieee_is_nan(this%spill_coefficient_of_discharge%value)) status = NML_ERR_NOT_SET
    case ("spill_coefficient_of_discharge%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_coefficient_of_discharge'"
        return
      end if
    case ("spill_coefficient_of_discharge%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_coefficient_of_discharge'"
        return
      end if
    case ("spill_coefficient_of_discharge%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_coefficient_of_discharge'"
        return
      end if
    case ("spill_coefficient_of_discharge")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_coefficient_of_discharge'"
        return
      end if
      if (ieee_is_nan(this%spill_coefficient_of_discharge%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("spill_exponent%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_exponent'"
        return
      end if
      if (ieee_is_nan(this%spill_exponent%value)) status = NML_ERR_NOT_SET
    case ("spill_exponent%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_exponent'"
        return
      end if
    case ("spill_exponent%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_exponent'"
        return
      end if
    case ("spill_exponent%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_exponent'"
        return
      end if
    case ("spill_exponent")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'spill_exponent'"
        return
      end if
      if (ieee_is_nan(this%spill_exponent%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("b%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'b'"
        return
      end if
      if (ieee_is_nan(this%b%value)) status = NML_ERR_NOT_SET
    case ("b%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'b'"
        return
      end if
    case ("b%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'b'"
        return
      end if
    case ("b%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'b'"
        return
      end if
    case ("b")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'b'"
        return
      end if
      if (ieee_is_nan(this%b%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_constant%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_constant%value)) status = NML_ERR_NOT_SET
    case ("ptf_ks_constant%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
    case ("ptf_ks_constant%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
    case ("ptf_ks_constant%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
    case ("ptf_ks_constant")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_constant'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_constant%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_sand%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_sand%value)) status = NML_ERR_NOT_SET
    case ("ptf_ks_sand%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
    case ("ptf_ks_sand%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
    case ("ptf_ks_sand%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
    case ("ptf_ks_sand")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_sand'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_sand%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("ptf_ks_clay%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_clay%value)) status = NML_ERR_NOT_SET
    case ("ptf_ks_clay%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
    case ("ptf_ks_clay%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
    case ("ptf_ks_clay%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
    case ("ptf_ks_clay")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'ptf_ks_clay'"
        return
      end if
      if (ieee_is_nan(this%ptf_ks_clay%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("baseflow_recession%value")
      if (.not. allocated(this%baseflow_recession)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%baseflow_recession), ubound(this%baseflow_recession), &
          "baseflow_recession", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%baseflow_recession(idx(1))%value)) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%baseflow_recession%value))) status = NML_ERR_NOT_SET
      end if
    case ("baseflow_recession%optimize")
      if (.not. allocated(this%baseflow_recession)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%baseflow_recession), ubound(this%baseflow_recession), &
          "baseflow_recession", errmsg)
        if (status /= NML_OK) return
      end if
    case ("baseflow_recession%min")
      if (.not. allocated(this%baseflow_recession)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%baseflow_recession), ubound(this%baseflow_recession), &
          "baseflow_recession", errmsg)
        if (status /= NML_OK) return
      end if
    case ("baseflow_recession%max")
      if (.not. allocated(this%baseflow_recession)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%baseflow_recession), ubound(this%baseflow_recession), &
          "baseflow_recession", errmsg)
        if (status /= NML_OK) return
      end if
    case ("baseflow_recession")
      if (.not. allocated(this%baseflow_recession)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%baseflow_recession), ubound(this%baseflow_recession), &
          "baseflow_recession", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%baseflow_recession(idx(1))%value)) then
          status = NML_ERR_NOT_SET
        end if
      else
        if (all(ieee_is_nan(this%baseflow_recession%value))) then
          status = NML_ERR_NOT_SET
        else if (any(ieee_is_nan(this%baseflow_recession%value))) then
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
    istat = this%is_set("b", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: b"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_ks_constant", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_ks_constant"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_ks_sand", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_ks_sand"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("ptf_ks_clay", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: ptf_ks_clay"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("baseflow_recession", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: baseflow_recession"
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
