!> \file nml_parameter_interflow_1.f90
!> \copydoc nml_interflow_1

!> \brief Interflow - Case 1
!> \details Parameters for parallel interflow reservoirs.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_interflow_1
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
    to_lower
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp
  use mo_parameter_types, only: parameter_t

  implicit none

  !> \class nml_interflow_1_t
  !> \brief Interflow - Case 1
  !> \details Parameters for parallel interflow reservoirs.
  type, public :: nml_interflow_1_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    type(parameter_t) :: storage_capacity_factor !< Interflow storage-capacity factor
    type(parameter_t) :: recession_slope !< Slope multiplier for interflow recession
    type(parameter_t) :: fast_recession_forest !< Forest multiplier for fast interflow recession
    type(parameter_t) :: slow_recession_ks !< Saturated-conductivity multiplier for slow interflow recession
    type(parameter_t) :: slow_recession_exponent !< Slow interflow exponent
  contains
    procedure :: init => nml_interflow_1_init
    procedure :: init_type => nml_interflow_1_init_type
    procedure :: from_file => nml_interflow_1_from_file
    procedure :: set => nml_interflow_1_set
    procedure :: is_set => nml_interflow_1_is_set
    procedure :: is_valid => nml_interflow_1_is_valid
  end type nml_interflow_1_t

contains

  !> \brief Initialize defaults and sentinels for interflow_1
  integer function nml_interflow_1_init(this, errmsg) result(status)
    class(nml_interflow_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! derived values
    status = this%init_type( &
      storage_capacity_factor=this%storage_capacity_factor, &
      recession_slope=this%recession_slope, &
      fast_recession_forest=this%fast_recession_forest, &
      slow_recession_ks=this%slow_recession_ks, &
      slow_recession_exponent=this%slow_recession_exponent, &
      errmsg=errmsg)
    if (status /= NML_OK) return
  end function nml_interflow_1_init

  !> \brief Initialize derived values with their field-specific defaults
  integer function nml_interflow_1_init_type(this, &
    storage_capacity_factor, &
    recession_slope, &
    fast_recession_forest, &
    slow_recession_ks, &
    slow_recession_exponent, &
    errmsg) result(status)
    class(nml_interflow_1_t), intent(in) :: this !< parent namelist instance
    type(parameter_t), intent(inout), optional :: storage_capacity_factor !< Interflow storage-capacity factor
    type(parameter_t), intent(inout), optional :: recession_slope !< Slope multiplier for interflow recession
    type(parameter_t), intent(inout), optional :: fast_recession_forest !< Forest multiplier for fast interflow recession
    type(parameter_t), intent(inout), optional :: slow_recession_ks !< Saturated-conductivity multiplier for slow interflow recession
    type(parameter_t), intent(inout), optional :: slow_recession_exponent !< Slow interflow exponent
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (present(storage_capacity_factor)) then
      storage_capacity_factor%value = ieee_value(storage_capacity_factor%value, ieee_quiet_nan) ! sentinel for derived component value
      storage_capacity_factor%optimize = .false.
      storage_capacity_factor%min = ieee_value(storage_capacity_factor%min, ieee_quiet_nan) ! sentinel for derived component min
      storage_capacity_factor%max = ieee_value(storage_capacity_factor%max, ieee_quiet_nan) ! sentinel for derived component max
      storage_capacity_factor%min = 75.0_dp
      storage_capacity_factor%max = 200.0_dp
    end if
    if (present(recession_slope)) then
      recession_slope%value = ieee_value(recession_slope%value, ieee_quiet_nan) ! sentinel for derived component value
      recession_slope%optimize = .false.
      recession_slope%min = ieee_value(recession_slope%min, ieee_quiet_nan) ! sentinel for derived component min
      recession_slope%max = ieee_value(recession_slope%max, ieee_quiet_nan) ! sentinel for derived component max
      recession_slope%min = 0.0_dp
      recession_slope%max = 10.0_dp
    end if
    if (present(fast_recession_forest)) then
      fast_recession_forest%value = ieee_value(fast_recession_forest%value, ieee_quiet_nan) ! sentinel for derived component value
      fast_recession_forest%optimize = .false.
      fast_recession_forest%min = ieee_value(fast_recession_forest%min, ieee_quiet_nan) ! sentinel for derived component min
      fast_recession_forest%max = ieee_value(fast_recession_forest%max, ieee_quiet_nan) ! sentinel for derived component max
      fast_recession_forest%min = 1.0_dp
      fast_recession_forest%max = 3.0_dp
    end if
    if (present(slow_recession_ks)) then
      slow_recession_ks%value = ieee_value(slow_recession_ks%value, ieee_quiet_nan) ! sentinel for derived component value
      slow_recession_ks%optimize = .false.
      slow_recession_ks%min = ieee_value(slow_recession_ks%min, ieee_quiet_nan) ! sentinel for derived component min
      slow_recession_ks%max = ieee_value(slow_recession_ks%max, ieee_quiet_nan) ! sentinel for derived component max
      slow_recession_ks%min = 1.0_dp
      slow_recession_ks%max = 30.0_dp
    end if
    if (present(slow_recession_exponent)) then
      slow_recession_exponent%value = ieee_value(slow_recession_exponent%value, ieee_quiet_nan) ! sentinel for derived component value
      slow_recession_exponent%optimize = .false.
      slow_recession_exponent%min = ieee_value(slow_recession_exponent%min, ieee_quiet_nan) ! sentinel for derived component min
      slow_recession_exponent%max = ieee_value(slow_recession_exponent%max, ieee_quiet_nan) ! sentinel for derived component max
      slow_recession_exponent%min = 0.05_dp
      slow_recession_exponent%max = 0.3_dp
    end if
  end function nml_interflow_1_init_type


  !> \brief Read interflow_1 namelist from file
  integer function nml_interflow_1_from_file(this, file, errmsg) result(status)
    class(nml_interflow_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    type(parameter_t) :: storage_capacity_factor
    type(parameter_t) :: recession_slope
    type(parameter_t) :: fast_recession_forest
    type(parameter_t) :: slow_recession_ks
    type(parameter_t) :: slow_recession_exponent
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /interflow_1/ &
      storage_capacity_factor, &
      recession_slope, &
      fast_recession_forest, &
      slow_recession_ks, &
      slow_recession_exponent

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    storage_capacity_factor = this%storage_capacity_factor
    recession_slope = this%recession_slope
    fast_recession_forest = this%fast_recession_forest
    slow_recession_ks = this%slow_recession_ks
    slow_recession_exponent = this%slow_recession_exponent

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("interflow_1", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=interflow_1, iostat=iostat, iomsg=iomsg)
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
    this%storage_capacity_factor = storage_capacity_factor
    this%recession_slope = recession_slope
    this%fast_recession_forest = fast_recession_forest
    this%slow_recession_ks = slow_recession_ks
    this%slow_recession_exponent = slow_recession_exponent

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_interflow_1_from_file

  !> \brief Set interflow_1 values
  integer function nml_interflow_1_set(this, &
    storage_capacity_factor, &
    recession_slope, &
    fast_recession_forest, &
    slow_recession_ks, &
    slow_recession_exponent, &
    errmsg) result(status)

    class(nml_interflow_1_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    type(parameter_t), intent(in) :: storage_capacity_factor !< Interflow storage-capacity factor
    type(parameter_t), intent(in) :: recession_slope !< Slope multiplier for interflow recession
    type(parameter_t), intent(in) :: fast_recession_forest !< Forest multiplier for fast interflow recession
    type(parameter_t), intent(in) :: slow_recession_ks !< Saturated-conductivity multiplier for slow interflow recession
    type(parameter_t), intent(in) :: slow_recession_exponent !< Slow interflow exponent

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    this%storage_capacity_factor = storage_capacity_factor
    this%recession_slope = recession_slope
    this%fast_recession_forest = fast_recession_forest
    this%slow_recession_ks = slow_recession_ks
    this%slow_recession_exponent = slow_recession_exponent

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_interflow_1_set

  !> \brief Check whether a namelist value was set
  integer function nml_interflow_1_is_set(this, name, idx, errmsg) result(status)
    class(nml_interflow_1_t), intent(in) :: this !< namelist instance
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
    case ("storage_capacity_factor%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'storage_capacity_factor'"
        return
      end if
      if (ieee_is_nan(this%storage_capacity_factor%value)) status = NML_ERR_NOT_SET
    case ("storage_capacity_factor%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'storage_capacity_factor'"
        return
      end if
    case ("storage_capacity_factor%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'storage_capacity_factor'"
        return
      end if
    case ("storage_capacity_factor%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'storage_capacity_factor'"
        return
      end if
    case ("storage_capacity_factor")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'storage_capacity_factor'"
        return
      end if
      if (ieee_is_nan(this%storage_capacity_factor%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("recession_slope%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recession_slope'"
        return
      end if
      if (ieee_is_nan(this%recession_slope%value)) status = NML_ERR_NOT_SET
    case ("recession_slope%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recession_slope'"
        return
      end if
    case ("recession_slope%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recession_slope'"
        return
      end if
    case ("recession_slope%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recession_slope'"
        return
      end if
    case ("recession_slope")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'recession_slope'"
        return
      end if
      if (ieee_is_nan(this%recession_slope%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("fast_recession_forest%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'fast_recession_forest'"
        return
      end if
      if (ieee_is_nan(this%fast_recession_forest%value)) status = NML_ERR_NOT_SET
    case ("fast_recession_forest%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'fast_recession_forest'"
        return
      end if
    case ("fast_recession_forest%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'fast_recession_forest'"
        return
      end if
    case ("fast_recession_forest%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'fast_recession_forest'"
        return
      end if
    case ("fast_recession_forest")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'fast_recession_forest'"
        return
      end if
      if (ieee_is_nan(this%fast_recession_forest%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("slow_recession_ks%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_ks'"
        return
      end if
      if (ieee_is_nan(this%slow_recession_ks%value)) status = NML_ERR_NOT_SET
    case ("slow_recession_ks%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_ks'"
        return
      end if
    case ("slow_recession_ks%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_ks'"
        return
      end if
    case ("slow_recession_ks%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_ks'"
        return
      end if
    case ("slow_recession_ks")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_ks'"
        return
      end if
      if (ieee_is_nan(this%slow_recession_ks%value)) then
        status = NML_ERR_NOT_SET
      end if
    case ("slow_recession_exponent%value")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_exponent'"
        return
      end if
      if (ieee_is_nan(this%slow_recession_exponent%value)) status = NML_ERR_NOT_SET
    case ("slow_recession_exponent%optimize")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_exponent'"
        return
      end if
    case ("slow_recession_exponent%min")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_exponent'"
        return
      end if
    case ("slow_recession_exponent%max")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_exponent'"
        return
      end if
    case ("slow_recession_exponent")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'slow_recession_exponent'"
        return
      end if
      if (ieee_is_nan(this%slow_recession_exponent%value)) then
        status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_interflow_1_is_set

  !> \brief Validate required values and constraints
  integer function nml_interflow_1_is_valid(this, errmsg) result(status)
    class(nml_interflow_1_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("storage_capacity_factor", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: storage_capacity_factor"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("recession_slope", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: recession_slope"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("fast_recession_forest", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: fast_recession_forest"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("slow_recession_ks", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: slow_recession_ks"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
    istat = this%is_set("slow_recession_exponent", errmsg=errmsg)
    if (istat == NML_ERR_NOT_SET) then
      status = NML_ERR_REQUIRED
      if (present(errmsg)) then
        if (len_trim(errmsg) == 0) then
          errmsg = "field not set: slow_recession_exponent"
        end if
        errmsg = "required " // trim(errmsg)
      end if
      return
    end if
    if (istat /= NML_OK) then
      status = istat
      return
    end if
  end function nml_interflow_1_is_valid

end module nml_interflow_1
