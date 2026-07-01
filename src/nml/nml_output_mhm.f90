!> \file nml_output_mhm.f90
!> \copydoc nml_output_mhm

!> \brief mHM output configuration
!> \details Output configuration for mHM.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_output_mhm
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
  integer(i4), parameter, public :: output_deflate_level__default = 6_i4
  logical, parameter, public :: output_double_precision__default = .false.
  integer(i4), parameter, public :: output_time_reference__default = 2_i4
  integer(i4), parameter, public :: output_frequency__default = -2_i4
  logical, parameter, public :: out_interception__default = .false.
  logical, parameter, public :: out_snowpack__default = .false.
  logical, parameter, public :: out_swc__default = .false.
  logical, parameter, public :: out_sm__default = .false.
  logical, parameter, public :: out_sm_all__default = .false.
  logical, parameter, public :: out_sealedstw__default = .false.
  logical, parameter, public :: out_unsatstw__default = .false.
  logical, parameter, public :: out_satstw__default = .false.
  logical, parameter, public :: out_pet__default = .false.
  logical, parameter, public :: out_aet_all__default = .false.
  logical, parameter, public :: out_q__default = .false.
  logical, parameter, public :: out_qd__default = .false.
  logical, parameter, public :: out_qif__default = .false.
  logical, parameter, public :: out_qis__default = .false.
  logical, parameter, public :: out_qb__default = .false.
  logical, parameter, public :: out_recharge__default = .false.
  logical, parameter, public :: out_soil_infil__default = .false.
  logical, parameter, public :: out_neutrons__default = .false.
  logical, parameter, public :: out_aet_layer__default = .false.
  logical, parameter, public :: out_preeffect__default = .false.
  logical, parameter, public :: out_qsm__default = .false.

  ! enum values
  integer(i4), parameter, public :: output_time_reference__enum_values(3) = [0_i4, 1_i4, 2_i4]

  ! bounds values
  integer(i4), parameter, public :: output_deflate_level__min = 0_i4
  integer(i4), parameter, public :: output_deflate_level__max = 9_i4
  integer(i4), parameter, public :: output_frequency__min = -3_i4

  !> \class nml_output_mhm_t
  !> \brief mHM output configuration
  !> \details Output configuration for mHM.
  type, public :: nml_output_mhm_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer(i4) :: output_deflate_level !< Output deflate level
    logical :: output_double_precision !< Output double precision
    integer(i4) :: output_time_reference !< Output time reference
    integer(i4) :: output_frequency !< Output time step
    logical :: out_interception !< Interception
    logical :: out_snowpack !< Snowpack
    logical :: out_swc !< Layered Soil Water Content
    logical :: out_sm !< Layered Volumetric Soil Moisture
    logical :: out_sm_all !< Mean Volumetric Soil Moisture
    logical :: out_sealedstw !< Reservoir of Sealed areas
    logical :: out_unsatstw !< Reservoir of Unsaturated areas
    logical :: out_satstw !< Reservoir of Saturated areas
    logical :: out_pet !< Potential Evapotranspiration
    logical :: out_aet_all !< Mean actual Evapotranspiration
    logical :: out_q !< Total Discharge
    logical :: out_qd !< Direct Runoff
    logical :: out_qif !< Fast Interflow
    logical :: out_qis !< Slow Interflow
    logical :: out_qb !< Baseflow
    logical :: out_recharge !< Groundwater Recharge
    logical :: out_soil_infil !< Soil Infiltration
    logical :: out_neutrons !< Neutrons
    logical :: out_aet_layer !< Actual Evapotranspiration from Soil Layers
    logical :: out_preeffect !< Effective Precipitation
    logical :: out_qsm !< Snow Melt
  contains
    procedure :: init => nml_output_mhm_init
    procedure :: from_file => nml_output_mhm_from_file
    procedure :: set => nml_output_mhm_set
    procedure :: is_set => nml_output_mhm_is_set
    procedure :: is_valid => nml_output_mhm_is_valid
  end type nml_output_mhm_t

contains

  !> \brief Check whether a value is part of an enum
  elemental logical function output_time_reference__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == output_time_reference__enum_values)
  end function output_time_reference__in_enum

  !> \brief Check whether a value is within bounds
  elemental logical function output_deflate_level__in_bounds(val, allow_missing) result(in_bounds)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_bounds = .true.
          return
        end if
      end if
    end if

    in_bounds = .true.
    if (val < output_deflate_level__min) in_bounds = .false.
    if (val > output_deflate_level__max) in_bounds = .false.
  end function output_deflate_level__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function output_frequency__in_bounds(val, allow_missing) result(in_bounds)
    integer(i4), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (val == -huge(val)) then
          in_bounds = .true.
          return
        end if
      end if
    end if

    in_bounds = .true.
    if (val < output_frequency__min) in_bounds = .false.
  end function output_frequency__in_bounds

  !> \brief Initialize defaults and sentinels for output_mhm
  integer function nml_output_mhm_init(this, errmsg) result(status)
    class(nml_output_mhm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! default values
    this%output_deflate_level = output_deflate_level__default
    this%output_double_precision = output_double_precision__default ! bool values always need a default
    this%output_time_reference = output_time_reference__default
    this%output_frequency = output_frequency__default
    this%out_interception = out_interception__default ! bool values always need a default
    this%out_snowpack = out_snowpack__default ! bool values always need a default
    this%out_swc = out_swc__default ! bool values always need a default
    this%out_sm = out_sm__default ! bool values always need a default
    this%out_sm_all = out_sm_all__default ! bool values always need a default
    this%out_sealedstw = out_sealedstw__default ! bool values always need a default
    this%out_unsatstw = out_unsatstw__default ! bool values always need a default
    this%out_satstw = out_satstw__default ! bool values always need a default
    this%out_pet = out_pet__default ! bool values always need a default
    this%out_aet_all = out_aet_all__default ! bool values always need a default
    this%out_q = out_q__default ! bool values always need a default
    this%out_qd = out_qd__default ! bool values always need a default
    this%out_qif = out_qif__default ! bool values always need a default
    this%out_qis = out_qis__default ! bool values always need a default
    this%out_qb = out_qb__default ! bool values always need a default
    this%out_recharge = out_recharge__default ! bool values always need a default
    this%out_soil_infil = out_soil_infil__default ! bool values always need a default
    this%out_neutrons = out_neutrons__default ! bool values always need a default
    this%out_aet_layer = out_aet_layer__default ! bool values always need a default
    this%out_preeffect = out_preeffect__default ! bool values always need a default
    this%out_qsm = out_qsm__default ! bool values always need a default
  end function nml_output_mhm_init


  !> \brief Read output_mhm namelist from file
  integer function nml_output_mhm_from_file(this, file, errmsg) result(status)
    class(nml_output_mhm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    integer(i4) :: output_deflate_level
    logical :: output_double_precision
    integer(i4) :: output_time_reference
    integer(i4) :: output_frequency
    logical :: out_interception
    logical :: out_snowpack
    logical :: out_swc
    logical :: out_sm
    logical :: out_sm_all
    logical :: out_sealedstw
    logical :: out_unsatstw
    logical :: out_satstw
    logical :: out_pet
    logical :: out_aet_all
    logical :: out_q
    logical :: out_qd
    logical :: out_qif
    logical :: out_qis
    logical :: out_qb
    logical :: out_recharge
    logical :: out_soil_infil
    logical :: out_neutrons
    logical :: out_aet_layer
    logical :: out_preeffect
    logical :: out_qsm
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /output_mhm/ &
      output_deflate_level, &
      output_double_precision, &
      output_time_reference, &
      output_frequency, &
      out_interception, &
      out_snowpack, &
      out_swc, &
      out_sm, &
      out_sm_all, &
      out_sealedstw, &
      out_unsatstw, &
      out_satstw, &
      out_pet, &
      out_aet_all, &
      out_q, &
      out_qd, &
      out_qif, &
      out_qis, &
      out_qb, &
      out_recharge, &
      out_soil_infil, &
      out_neutrons, &
      out_aet_layer, &
      out_preeffect, &
      out_qsm

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    output_deflate_level = this%output_deflate_level
    output_double_precision = this%output_double_precision
    output_time_reference = this%output_time_reference
    output_frequency = this%output_frequency
    out_interception = this%out_interception
    out_snowpack = this%out_snowpack
    out_swc = this%out_swc
    out_sm = this%out_sm
    out_sm_all = this%out_sm_all
    out_sealedstw = this%out_sealedstw
    out_unsatstw = this%out_unsatstw
    out_satstw = this%out_satstw
    out_pet = this%out_pet
    out_aet_all = this%out_aet_all
    out_q = this%out_q
    out_qd = this%out_qd
    out_qif = this%out_qif
    out_qis = this%out_qis
    out_qb = this%out_qb
    out_recharge = this%out_recharge
    out_soil_infil = this%out_soil_infil
    out_neutrons = this%out_neutrons
    out_aet_layer = this%out_aet_layer
    out_preeffect = this%out_preeffect
    out_qsm = this%out_qsm

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("output_mhm", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=output_mhm, iostat=iostat, iomsg=iomsg)
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
    this%output_deflate_level = output_deflate_level
    this%output_double_precision = output_double_precision
    this%output_time_reference = output_time_reference
    this%output_frequency = output_frequency
    this%out_interception = out_interception
    this%out_snowpack = out_snowpack
    this%out_swc = out_swc
    this%out_sm = out_sm
    this%out_sm_all = out_sm_all
    this%out_sealedstw = out_sealedstw
    this%out_unsatstw = out_unsatstw
    this%out_satstw = out_satstw
    this%out_pet = out_pet
    this%out_aet_all = out_aet_all
    this%out_q = out_q
    this%out_qd = out_qd
    this%out_qif = out_qif
    this%out_qis = out_qis
    this%out_qb = out_qb
    this%out_recharge = out_recharge
    this%out_soil_infil = out_soil_infil
    this%out_neutrons = out_neutrons
    this%out_aet_layer = out_aet_layer
    this%out_preeffect = out_preeffect
    this%out_qsm = out_qsm

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_output_mhm_from_file

  !> \brief Set output_mhm values
  integer function nml_output_mhm_set(this, &
    output_deflate_level, &
    output_double_precision, &
    output_time_reference, &
    output_frequency, &
    out_interception, &
    out_snowpack, &
    out_swc, &
    out_sm, &
    out_sm_all, &
    out_sealedstw, &
    out_unsatstw, &
    out_satstw, &
    out_pet, &
    out_aet_all, &
    out_q, &
    out_qd, &
    out_qif, &
    out_qis, &
    out_qb, &
    out_recharge, &
    out_soil_infil, &
    out_neutrons, &
    out_aet_layer, &
    out_preeffect, &
    out_qsm, &
    errmsg) result(status)

    class(nml_output_mhm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer(i4), intent(in), optional :: output_deflate_level !< Output deflate level
    logical, intent(in), optional :: output_double_precision !< Output double precision
    integer(i4), intent(in), optional :: output_time_reference !< Output time reference
    integer(i4), intent(in), optional :: output_frequency !< Output time step
    logical, intent(in), optional :: out_interception !< Interception
    logical, intent(in), optional :: out_snowpack !< Snowpack
    logical, intent(in), optional :: out_swc !< Layered Soil Water Content
    logical, intent(in), optional :: out_sm !< Layered Volumetric Soil Moisture
    logical, intent(in), optional :: out_sm_all !< Mean Volumetric Soil Moisture
    logical, intent(in), optional :: out_sealedstw !< Reservoir of Sealed areas
    logical, intent(in), optional :: out_unsatstw !< Reservoir of Unsaturated areas
    logical, intent(in), optional :: out_satstw !< Reservoir of Saturated areas
    logical, intent(in), optional :: out_pet !< Potential Evapotranspiration
    logical, intent(in), optional :: out_aet_all !< Mean actual Evapotranspiration
    logical, intent(in), optional :: out_q !< Total Discharge
    logical, intent(in), optional :: out_qd !< Direct Runoff
    logical, intent(in), optional :: out_qif !< Fast Interflow
    logical, intent(in), optional :: out_qis !< Slow Interflow
    logical, intent(in), optional :: out_qb !< Baseflow
    logical, intent(in), optional :: out_recharge !< Groundwater Recharge
    logical, intent(in), optional :: out_soil_infil !< Soil Infiltration
    logical, intent(in), optional :: out_neutrons !< Neutrons
    logical, intent(in), optional :: out_aet_layer !< Actual Evapotranspiration from Soil Layers
    logical, intent(in), optional :: out_preeffect !< Effective Precipitation
    logical, intent(in), optional :: out_qsm !< Snow Melt

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(output_deflate_level)) this%output_deflate_level = output_deflate_level
    if (present(output_double_precision)) this%output_double_precision = output_double_precision
    if (present(output_time_reference)) this%output_time_reference = output_time_reference
    if (present(output_frequency)) this%output_frequency = output_frequency
    if (present(out_interception)) this%out_interception = out_interception
    if (present(out_snowpack)) this%out_snowpack = out_snowpack
    if (present(out_swc)) this%out_swc = out_swc
    if (present(out_sm)) this%out_sm = out_sm
    if (present(out_sm_all)) this%out_sm_all = out_sm_all
    if (present(out_sealedstw)) this%out_sealedstw = out_sealedstw
    if (present(out_unsatstw)) this%out_unsatstw = out_unsatstw
    if (present(out_satstw)) this%out_satstw = out_satstw
    if (present(out_pet)) this%out_pet = out_pet
    if (present(out_aet_all)) this%out_aet_all = out_aet_all
    if (present(out_q)) this%out_q = out_q
    if (present(out_qd)) this%out_qd = out_qd
    if (present(out_qif)) this%out_qif = out_qif
    if (present(out_qis)) this%out_qis = out_qis
    if (present(out_qb)) this%out_qb = out_qb
    if (present(out_recharge)) this%out_recharge = out_recharge
    if (present(out_soil_infil)) this%out_soil_infil = out_soil_infil
    if (present(out_neutrons)) this%out_neutrons = out_neutrons
    if (present(out_aet_layer)) this%out_aet_layer = out_aet_layer
    if (present(out_preeffect)) this%out_preeffect = out_preeffect
    if (present(out_qsm)) this%out_qsm = out_qsm

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_output_mhm_set

  !> \brief Check whether a namelist value was set
  integer function nml_output_mhm_is_set(this, name, idx, errmsg) result(status)
    class(nml_output_mhm_t), intent(in) :: this !< namelist instance
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
    case ("output_deflate_level")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'output_deflate_level'"
        return
      end if
    case ("output_double_precision")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'output_double_precision'"
        return
      end if
    case ("output_time_reference")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'output_time_reference'"
        return
      end if
    case ("output_frequency")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'output_frequency'"
        return
      end if
    case ("out_interception")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_interception'"
        return
      end if
    case ("out_snowpack")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_snowpack'"
        return
      end if
    case ("out_swc")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_SWC'"
        return
      end if
    case ("out_sm")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_SM'"
        return
      end if
    case ("out_sm_all")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_SM_all'"
        return
      end if
    case ("out_sealedstw")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_sealedSTW'"
        return
      end if
    case ("out_unsatstw")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_unsatSTW'"
        return
      end if
    case ("out_satstw")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_satSTW'"
        return
      end if
    case ("out_pet")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_PET'"
        return
      end if
    case ("out_aet_all")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_aET_all'"
        return
      end if
    case ("out_q")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_Q'"
        return
      end if
    case ("out_qd")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_QD'"
        return
      end if
    case ("out_qif")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_QIf'"
        return
      end if
    case ("out_qis")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_QIs'"
        return
      end if
    case ("out_qb")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_QB'"
        return
      end if
    case ("out_recharge")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_recharge'"
        return
      end if
    case ("out_soil_infil")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_soil_infil'"
        return
      end if
    case ("out_neutrons")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_neutrons'"
        return
      end if
    case ("out_aet_layer")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_aET_layer'"
        return
      end if
    case ("out_preeffect")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_preEffect'"
        return
      end if
    case ("out_qsm")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'out_Qsm'"
        return
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_output_mhm_is_set

  !> \brief Validate required values and constraints
  integer function nml_output_mhm_is_valid(this, errmsg) result(status)
    class(nml_output_mhm_t), intent(in) :: this !< namelist instance
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
    istat = this%is_set("output_time_reference", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. output_time_reference__in_enum(this%output_time_reference)) then
        status = NML_ERR_ENUM
        if (present(errmsg)) errmsg = "enum constraint failed: output_time_reference"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    ! bounds constraints
    istat = this%is_set("output_deflate_level", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. output_deflate_level__in_bounds(this%output_deflate_level)) then
        status = NML_ERR_BOUNDS
        if (present(errmsg)) errmsg = "bounds constraint failed: output_deflate_level"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
    istat = this%is_set("output_frequency", errmsg=errmsg)
    if (istat == NML_OK) then
      if (.not. output_frequency__in_bounds(this%output_frequency)) then
        status = NML_ERR_BOUNDS
        if (present(errmsg)) errmsg = "bounds constraint failed: output_frequency"
        return
      end if
    else if (istat /= NML_ERR_NOT_SET) then
      status = istat
      return
    end if
  end function nml_output_mhm_is_valid

end module nml_output_mhm
