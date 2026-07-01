!> \file nml_config_processes.f90
!> \copydoc nml_config_processes

!> \brief Processes configuration
!> \details Configuration for process case selection in mHM.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_processes
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
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    i4

  implicit none

  ! default values
  integer(i4), parameter, public :: interception__default = 0_i4
  integer(i4), parameter, public :: snow__default = 0_i4
  integer(i4), parameter, public :: soil_moisture__default = 0_i4
  integer(i4), parameter, public :: direct_runoff__default = 0_i4
  integer(i4), parameter, public :: pet__default = 0_i4
  integer(i4), parameter, public :: interflow__default = 0_i4
  integer(i4), parameter, public :: percolation__default = 0_i4
  integer(i4), parameter, public :: baseflow__default = 0_i4
  integer(i4), parameter, public :: neutrons__default = 0_i4
  integer(i4), parameter, public :: routing__default = 0_i4
  integer(i4), parameter, public :: temperature_routing__default = 0_i4

  ! enum values
  integer(i4), parameter, public :: interception__enum_values(3) = [-1_i4, 0_i4, 1_i4]
  integer(i4), parameter, public :: snow__enum_values(3) = [-1_i4, 0_i4, 1_i4]
  integer(i4), parameter, public :: soil_moisture__enum_values(5) = [0_i4, 1_i4, 2_i4, 3_i4, 4_i4]
  integer(i4), parameter, public :: direct_runoff__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: pet__enum_values(6) = [-2_i4, -1_i4, 0_i4, 1_i4, 2_i4, 3_i4]
  integer(i4), parameter, public :: interflow__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: percolation__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: baseflow__enum_values(2) = [0_i4, 1_i4]
  integer(i4), parameter, public :: neutrons__enum_values(3) = [0_i4, 1_i4, 2_i4]
  integer(i4), parameter, public :: routing__enum_values(4) = [0_i4, 1_i4, 2_i4, 3_i4]
  integer(i4), parameter, public :: temperature_routing__enum_values(2) = [0_i4, 1_i4]

  !> \class nml_config_processes_t
  !> \brief Processes configuration
  !> \details Configuration for process case selection in mHM.
  type, public :: nml_config_processes_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer(i4) :: interception !< Interception process case
    integer(i4) :: snow !< Snow process case
    integer(i4) :: soil_moisture !< Soil moisture process case
    integer(i4) :: direct_runoff !< Direct runoff process case
    integer(i4) :: pet !< Potential evapotranspiration (PET) process case
    integer(i4) :: interflow !< Interflow process case
    integer(i4) :: percolation !< Percolation process case
    integer(i4) :: baseflow !< Baseflow process case
    integer(i4) :: neutrons !< Ground albedo of cosmic-ray neutrons process case
    integer(i4) :: routing !< Routing process case
    integer(i4) :: temperature_routing !< River temperature routing process case
  contains
    procedure :: init => nml_config_processes_init
    procedure :: from_file => nml_config_processes_from_file
    procedure :: set => nml_config_processes_set
    procedure :: is_set => nml_config_processes_is_set
    procedure :: is_valid => nml_config_processes_is_valid
  end type nml_config_processes_t

contains

  !> \brief Check whether a value is part of an enum
  elemental logical function interception__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == interception__enum_values)
  end function interception__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function snow__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == snow__enum_values)
  end function snow__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function soil_moisture__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == soil_moisture__enum_values)
  end function soil_moisture__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function direct_runoff__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == direct_runoff__enum_values)
  end function direct_runoff__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function pet__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == pet__enum_values)
  end function pet__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function interflow__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == interflow__enum_values)
  end function interflow__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function percolation__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == percolation__enum_values)
  end function percolation__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function baseflow__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == baseflow__enum_values)
  end function baseflow__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function neutrons__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == neutrons__enum_values)
  end function neutrons__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function routing__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == routing__enum_values)
  end function routing__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function temperature_routing__in_enum(val, allow_missing) result(in_enum)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_enum = .true.
          return
        end if
      end if
    end if
    in_enum = any(val == temperature_routing__enum_values)
  end function temperature_routing__in_enum

  !> \brief Initialize defaults and sentinels for config_processes
  integer function nml_config_processes_init(this, errmsg) result(status)
    class(nml_config_processes_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! default values
    this%interception = interception__default
    this%snow = snow__default
    this%soil_moisture = soil_moisture__default
    this%direct_runoff = direct_runoff__default
    this%pet = pet__default
    this%interflow = interflow__default
    this%percolation = percolation__default
    this%baseflow = baseflow__default
    this%neutrons = neutrons__default
    this%routing = routing__default
    this%temperature_routing = temperature_routing__default
  end function nml_config_processes_init


  !> \brief Read config_processes namelist from file
  integer function nml_config_processes_from_file(this, file, errmsg) result(status)
    class(nml_config_processes_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    integer(i4) :: interception
    integer(i4) :: snow
    integer(i4) :: soil_moisture
    integer(i4) :: direct_runoff
    integer(i4) :: pet
    integer(i4) :: interflow
    integer(i4) :: percolation
    integer(i4) :: baseflow
    integer(i4) :: neutrons
    integer(i4) :: routing
    integer(i4) :: temperature_routing
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_processes/ &
      interception, &
      snow, &
      soil_moisture, &
      direct_runoff, &
      pet, &
      interflow, &
      percolation, &
      baseflow, &
      neutrons, &
      routing, &
      temperature_routing

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    interception = this%interception
    snow = this%snow
    soil_moisture = this%soil_moisture
    direct_runoff = this%direct_runoff
    pet = this%pet
    interflow = this%interflow
    percolation = this%percolation
    baseflow = this%baseflow
    neutrons = this%neutrons
    routing = this%routing
    temperature_routing = this%temperature_routing

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_processes", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_processes, iostat=iostat, iomsg=iomsg)
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
    this%interception = interception
    this%snow = snow
    this%soil_moisture = soil_moisture
    this%direct_runoff = direct_runoff
    this%pet = pet
    this%interflow = interflow
    this%percolation = percolation
    this%baseflow = baseflow
    this%neutrons = neutrons
    this%routing = routing
    this%temperature_routing = temperature_routing

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_processes_from_file

  !> \brief Set config_processes values
  integer function nml_config_processes_set(this, &
    interception, &
    snow, &
    soil_moisture, &
    direct_runoff, &
    pet, &
    interflow, &
    percolation, &
    baseflow, &
    neutrons, &
    routing, &
    temperature_routing, &
    errmsg) result(status)

    class(nml_config_processes_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer(i4), intent(in), optional :: interception !< Interception process case
    integer(i4), intent(in), optional :: snow !< Snow process case
    integer(i4), intent(in), optional :: soil_moisture !< Soil moisture process case
    integer(i4), intent(in), optional :: direct_runoff !< Direct runoff process case
    integer(i4), intent(in), optional :: pet !< Potential evapotranspiration (PET) process case
    integer(i4), intent(in), optional :: interflow !< Interflow process case
    integer(i4), intent(in), optional :: percolation !< Percolation process case
    integer(i4), intent(in), optional :: baseflow !< Baseflow process case
    integer(i4), intent(in), optional :: neutrons !< Ground albedo of cosmic-ray neutrons process case
    integer(i4), intent(in), optional :: routing !< Routing process case
    integer(i4), intent(in), optional :: temperature_routing !< River temperature routing process case

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(interception)) this%interception = interception
    if (present(snow)) this%snow = snow
    if (present(soil_moisture)) this%soil_moisture = soil_moisture
    if (present(direct_runoff)) this%direct_runoff = direct_runoff
    if (present(pet)) this%pet = pet
    if (present(interflow)) this%interflow = interflow
    if (present(percolation)) this%percolation = percolation
    if (present(baseflow)) this%baseflow = baseflow
    if (present(neutrons)) this%neutrons = neutrons
    if (present(routing)) this%routing = routing
    if (present(temperature_routing)) this%temperature_routing = temperature_routing

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_processes_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_processes_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_processes_t), intent(in) :: this !< namelist instance
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
    case ("interception")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'interception'"
        return
      end if
    case ("snow")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'snow'"
        return
      end if
    case ("soil_moisture")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'soil_moisture'"
        return
      end if
    case ("direct_runoff")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'direct_runoff'"
        return
      end if
    case ("pet")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'pet'"
        return
      end if
    case ("interflow")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'interflow'"
        return
      end if
    case ("percolation")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'percolation'"
        return
      end if
    case ("baseflow")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'baseflow'"
        return
      end if
    case ("neutrons")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'neutrons'"
        return
      end if
    case ("routing")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'routing'"
        return
      end if
    case ("temperature_routing")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'temperature_routing'"
        return
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_processes_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_processes_is_valid(this, errmsg) result(status)
    class(nml_config_processes_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

    ! enum constraints
    istat = this%is_set("interception", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. interception__in_enum(this%interception)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: interception"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("snow", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. snow__in_enum(this%snow)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: snow"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("soil_moisture", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. soil_moisture__in_enum(this%soil_moisture)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: soil_moisture"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("direct_runoff", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. direct_runoff__in_enum(this%direct_runoff)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: direct_runoff"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("pet", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. pet__in_enum(this%pet)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: pet"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("interflow", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. interflow__in_enum(this%interflow)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: interflow"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("percolation", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. percolation__in_enum(this%percolation)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: percolation"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("baseflow", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. baseflow__in_enum(this%baseflow)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: baseflow"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("neutrons", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. neutrons__in_enum(this%neutrons)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: neutrons"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("routing", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. routing__in_enum(this%routing)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: routing"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("temperature_routing", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. temperature_routing__in_enum(this%temperature_routing)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: temperature_routing"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
  end function nml_config_processes_is_valid

end module nml_config_processes
