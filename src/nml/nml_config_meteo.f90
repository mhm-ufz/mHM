!> \file nml_config_meteo.f90
!> \copydoc nml_config_meteo

!> \brief Meteorological configuration
!> \details Configuration for meteorological input data handling in mHM.
!! Meteorological weights can be used to disaggregate daily data to hourly values.
!! Fixed-hour forcing support must match the global model step. A 24-hour model
!! step therefore requires daily or fixed 24-hour forcing; aggregation of hourly
!! forcing to daily model steps is not implemented.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_meteo
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
    buf
  use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  ! kind specifiers listed in the nml-tools configuration file
  use mo_kind, only: &
    dp

  implicit none

  ! default values
  logical, parameter, public :: read_meteo_weights__default = .false.
  character(len=buf), parameter, public :: pre_weights_var__default = "pre_weights"
  character(len=buf), parameter, public :: pet_weights_var__default = "pet_weights"
  character(len=buf), parameter, public :: temp_weights_var__default = "tavg_weights"
  character(len=buf), parameter, public :: ssrd_weights_var__default = "ssrd_weights"
  character(len=buf), parameter, public :: strd_weights_var__default = "strd_weights"
  logical, parameter, public :: share_frac__default = .true.

  !> \class nml_config_meteo_t
  !> \brief Meteorological configuration
  !> \details Configuration for meteorological input data handling in mHM.
  !! Meteorological weights can be used to disaggregate daily data to hourly values.
  !! Fixed-hour forcing support must match the global model step. A 24-hour model
  !! step therefore requires daily or fixed 24-hour forcing; aggregation of hourly
  !! forcing to daily model steps is not implemented.
  type, public :: nml_config_meteo_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    logical, allocatable, dimension(:) :: read_meteo_weights !< Read meteorological weights
    character(len=buf), allocatable, dimension(:) :: pre_weights_path !< Precipitation weights path
    character(len=buf), allocatable, dimension(:) :: pet_weights_path !< Potential evapotranspiration weights path
    character(len=buf), allocatable, dimension(:) :: temp_weights_path !< Surface downward shortwave radiation weights path
    character(len=buf), allocatable, dimension(:) :: ssrd_weights_path !< Surface downward shortwave radiation weights path
    character(len=buf), allocatable, dimension(:) :: strd_weights_path !< Surface downward longwave radiation weights path
    character(len=buf), allocatable, dimension(:) :: pre_weights_var !< Precipitation weights variable name
    character(len=buf), allocatable, dimension(:) :: pet_weights_var !< Potential evapotranspiration weights variable name
    character(len=buf), allocatable, dimension(:) :: temp_weights_var !< Average temperature weights variable name
    character(len=buf), allocatable, dimension(:) :: ssrd_weights_var !< Surface downward shortwave radiation weights variable name
    character(len=buf), allocatable, dimension(:) :: strd_weights_var !< Surface downward longwave radiation weights variable name
    real(dp), allocatable, dimension(:, :) :: frac_night_pre !< Fraction of nightly precipitation
    real(dp), allocatable, dimension(:, :) :: frac_night_pet !< Fraction of nightly potential evapotranspiration
    real(dp), allocatable, dimension(:, :) :: frac_night_temp !< Fraction of nightly temperature
    real(dp), allocatable, dimension(:, :) :: frac_night_ssrd !< Fraction of nightly surface downward shortwave radiation
    real(dp), allocatable, dimension(:, :) :: frac_night_strd !< Fraction of nightly surface downward longwave radiation
    logical :: share_frac !< Share fractions between domains
  contains
    procedure :: init => nml_config_meteo_init
    procedure :: set_dims => nml_config_meteo_set_dims
    procedure :: from_file => nml_config_meteo_from_file
    procedure :: set => nml_config_meteo_set
    procedure :: is_set => nml_config_meteo_is_set
    procedure :: is_valid => nml_config_meteo_is_valid
  end type nml_config_meteo_t

contains

  !> \brief Initialize defaults and sentinels for config_meteo
  integer function nml_config_meteo_init(this, errmsg) result(status)
    class(nml_config_meteo_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%read_meteo_weights)) deallocate(this%read_meteo_weights)
    allocate(this%read_meteo_weights(this%n_domains))
    if (allocated(this%pre_weights_path)) deallocate(this%pre_weights_path)
    allocate(character(len=buf) :: this%pre_weights_path(this%n_domains))
    if (allocated(this%pet_weights_path)) deallocate(this%pet_weights_path)
    allocate(character(len=buf) :: this%pet_weights_path(this%n_domains))
    if (allocated(this%temp_weights_path)) deallocate(this%temp_weights_path)
    allocate(character(len=buf) :: this%temp_weights_path(this%n_domains))
    if (allocated(this%ssrd_weights_path)) deallocate(this%ssrd_weights_path)
    allocate(character(len=buf) :: this%ssrd_weights_path(this%n_domains))
    if (allocated(this%strd_weights_path)) deallocate(this%strd_weights_path)
    allocate(character(len=buf) :: this%strd_weights_path(this%n_domains))
    if (allocated(this%pre_weights_var)) deallocate(this%pre_weights_var)
    allocate(character(len=buf) :: this%pre_weights_var(this%n_domains))
    if (allocated(this%pet_weights_var)) deallocate(this%pet_weights_var)
    allocate(character(len=buf) :: this%pet_weights_var(this%n_domains))
    if (allocated(this%temp_weights_var)) deallocate(this%temp_weights_var)
    allocate(character(len=buf) :: this%temp_weights_var(this%n_domains))
    if (allocated(this%ssrd_weights_var)) deallocate(this%ssrd_weights_var)
    allocate(character(len=buf) :: this%ssrd_weights_var(this%n_domains))
    if (allocated(this%strd_weights_var)) deallocate(this%strd_weights_var)
    allocate(character(len=buf) :: this%strd_weights_var(this%n_domains))
    if (allocated(this%frac_night_pre)) deallocate(this%frac_night_pre)
    allocate(this%frac_night_pre(12, this%n_domains))
    if (allocated(this%frac_night_pet)) deallocate(this%frac_night_pet)
    allocate(this%frac_night_pet(12, this%n_domains))
    if (allocated(this%frac_night_temp)) deallocate(this%frac_night_temp)
    allocate(this%frac_night_temp(12, this%n_domains))
    if (allocated(this%frac_night_ssrd)) deallocate(this%frac_night_ssrd)
    allocate(this%frac_night_ssrd(12, this%n_domains))
    if (allocated(this%frac_night_strd)) deallocate(this%frac_night_strd)
    allocate(this%frac_night_strd(12, this%n_domains))

    ! sentinel values for required/optional parameters
    this%pre_weights_path = achar(0) ! sentinel for optional string array
    this%pet_weights_path = achar(0) ! sentinel for optional string array
    this%temp_weights_path = achar(0) ! sentinel for optional string array
    this%ssrd_weights_path = achar(0) ! sentinel for optional string array
    this%strd_weights_path = achar(0) ! sentinel for optional string array
    this%frac_night_pre = ieee_value(this%frac_night_pre, ieee_quiet_nan) ! sentinel for optional real array
    this%frac_night_pet = ieee_value(this%frac_night_pet, ieee_quiet_nan) ! sentinel for optional real array
    this%frac_night_temp = ieee_value(this%frac_night_temp, ieee_quiet_nan) ! sentinel for optional real array
    this%frac_night_ssrd = ieee_value(this%frac_night_ssrd, ieee_quiet_nan) ! sentinel for optional real array
    this%frac_night_strd = ieee_value(this%frac_night_strd, ieee_quiet_nan) ! sentinel for optional real array
    ! default values
    this%read_meteo_weights = read_meteo_weights__default
    this%pre_weights_var = pre_weights_var__default
    this%pet_weights_var = pet_weights_var__default
    this%temp_weights_var = temp_weights_var__default
    this%ssrd_weights_var = ssrd_weights_var__default
    this%strd_weights_var = strd_weights_var__default
    this%share_frac = share_frac__default ! bool values always need a default
  end function nml_config_meteo_init

  !> \brief Reset runtime dimensions for config_meteo
  integer function nml_config_meteo_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_meteo_t), intent(inout) :: this !< namelist instance
    integer, intent(in), optional :: n_domains !< runtime dimension override for n_domains
    integer :: candidate__n_domains
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
    this%n_domains = candidate__n_domains

    ! deallocate runtime-sized fields; init/set/from_file allocate them again
    if (allocated(this%read_meteo_weights)) deallocate(this%read_meteo_weights)
    if (allocated(this%pre_weights_path)) deallocate(this%pre_weights_path)
    if (allocated(this%pet_weights_path)) deallocate(this%pet_weights_path)
    if (allocated(this%temp_weights_path)) deallocate(this%temp_weights_path)
    if (allocated(this%ssrd_weights_path)) deallocate(this%ssrd_weights_path)
    if (allocated(this%strd_weights_path)) deallocate(this%strd_weights_path)
    if (allocated(this%pre_weights_var)) deallocate(this%pre_weights_var)
    if (allocated(this%pet_weights_var)) deallocate(this%pet_weights_var)
    if (allocated(this%temp_weights_var)) deallocate(this%temp_weights_var)
    if (allocated(this%ssrd_weights_var)) deallocate(this%ssrd_weights_var)
    if (allocated(this%strd_weights_var)) deallocate(this%strd_weights_var)
    if (allocated(this%frac_night_pre)) deallocate(this%frac_night_pre)
    if (allocated(this%frac_night_pet)) deallocate(this%frac_night_pet)
    if (allocated(this%frac_night_temp)) deallocate(this%frac_night_temp)
    if (allocated(this%frac_night_ssrd)) deallocate(this%frac_night_ssrd)
    if (allocated(this%frac_night_strd)) deallocate(this%frac_night_strd)
    this%is_configured = .false.
  end function nml_config_meteo_set_dims


  !> \brief Read config_meteo namelist from file
  integer function nml_config_meteo_from_file(this, file, errmsg) result(status)
    class(nml_config_meteo_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    logical, allocatable, dimension(:) :: read_meteo_weights
    character(len=buf), allocatable, dimension(:) :: pre_weights_path
    character(len=buf), allocatable, dimension(:) :: pet_weights_path
    character(len=buf), allocatable, dimension(:) :: temp_weights_path
    character(len=buf), allocatable, dimension(:) :: ssrd_weights_path
    character(len=buf), allocatable, dimension(:) :: strd_weights_path
    character(len=buf), allocatable, dimension(:) :: pre_weights_var
    character(len=buf), allocatable, dimension(:) :: pet_weights_var
    character(len=buf), allocatable, dimension(:) :: temp_weights_var
    character(len=buf), allocatable, dimension(:) :: ssrd_weights_var
    character(len=buf), allocatable, dimension(:) :: strd_weights_var
    real(dp), allocatable, dimension(:, :) :: frac_night_pre
    real(dp), allocatable, dimension(:, :) :: frac_night_pet
    real(dp), allocatable, dimension(:, :) :: frac_night_temp
    real(dp), allocatable, dimension(:, :) :: frac_night_ssrd
    real(dp), allocatable, dimension(:, :) :: frac_night_strd
    logical :: share_frac
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_meteo/ &
      read_meteo_weights, &
      pre_weights_path, &
      pet_weights_path, &
      temp_weights_path, &
      ssrd_weights_path, &
      strd_weights_path, &
      pre_weights_var, &
      pet_weights_var, &
      temp_weights_var, &
      ssrd_weights_var, &
      strd_weights_var, &
      frac_night_pre, &
      frac_night_pet, &
      frac_night_temp, &
      frac_night_ssrd, &
      frac_night_strd, &
      share_frac

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(read_meteo_weights)) deallocate(read_meteo_weights)
    allocate(read_meteo_weights(this%n_domains))
    if (allocated(pre_weights_path)) deallocate(pre_weights_path)
    allocate(character(len=buf) :: pre_weights_path(this%n_domains))
    if (allocated(pet_weights_path)) deallocate(pet_weights_path)
    allocate(character(len=buf) :: pet_weights_path(this%n_domains))
    if (allocated(temp_weights_path)) deallocate(temp_weights_path)
    allocate(character(len=buf) :: temp_weights_path(this%n_domains))
    if (allocated(ssrd_weights_path)) deallocate(ssrd_weights_path)
    allocate(character(len=buf) :: ssrd_weights_path(this%n_domains))
    if (allocated(strd_weights_path)) deallocate(strd_weights_path)
    allocate(character(len=buf) :: strd_weights_path(this%n_domains))
    if (allocated(pre_weights_var)) deallocate(pre_weights_var)
    allocate(character(len=buf) :: pre_weights_var(this%n_domains))
    if (allocated(pet_weights_var)) deallocate(pet_weights_var)
    allocate(character(len=buf) :: pet_weights_var(this%n_domains))
    if (allocated(temp_weights_var)) deallocate(temp_weights_var)
    allocate(character(len=buf) :: temp_weights_var(this%n_domains))
    if (allocated(ssrd_weights_var)) deallocate(ssrd_weights_var)
    allocate(character(len=buf) :: ssrd_weights_var(this%n_domains))
    if (allocated(strd_weights_var)) deallocate(strd_weights_var)
    allocate(character(len=buf) :: strd_weights_var(this%n_domains))
    if (allocated(frac_night_pre)) deallocate(frac_night_pre)
    allocate(frac_night_pre(12, this%n_domains))
    if (allocated(frac_night_pet)) deallocate(frac_night_pet)
    allocate(frac_night_pet(12, this%n_domains))
    if (allocated(frac_night_temp)) deallocate(frac_night_temp)
    allocate(frac_night_temp(12, this%n_domains))
    if (allocated(frac_night_ssrd)) deallocate(frac_night_ssrd)
    allocate(frac_night_ssrd(12, this%n_domains))
    if (allocated(frac_night_strd)) deallocate(frac_night_strd)
    allocate(frac_night_strd(12, this%n_domains))
    read_meteo_weights = this%read_meteo_weights
    pre_weights_path = this%pre_weights_path
    pet_weights_path = this%pet_weights_path
    temp_weights_path = this%temp_weights_path
    ssrd_weights_path = this%ssrd_weights_path
    strd_weights_path = this%strd_weights_path
    pre_weights_var = this%pre_weights_var
    pet_weights_var = this%pet_weights_var
    temp_weights_var = this%temp_weights_var
    ssrd_weights_var = this%ssrd_weights_var
    strd_weights_var = this%strd_weights_var
    frac_night_pre = this%frac_night_pre
    frac_night_pet = this%frac_night_pet
    frac_night_temp = this%frac_night_temp
    frac_night_ssrd = this%frac_night_ssrd
    frac_night_strd = this%frac_night_strd
    share_frac = this%share_frac

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_meteo", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_meteo, iostat=iostat, iomsg=iomsg)
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
    this%read_meteo_weights = read_meteo_weights
    this%pre_weights_path = pre_weights_path
    this%pet_weights_path = pet_weights_path
    this%temp_weights_path = temp_weights_path
    this%ssrd_weights_path = ssrd_weights_path
    this%strd_weights_path = strd_weights_path
    this%pre_weights_var = pre_weights_var
    this%pet_weights_var = pet_weights_var
    this%temp_weights_var = temp_weights_var
    this%ssrd_weights_var = ssrd_weights_var
    this%strd_weights_var = strd_weights_var
    this%frac_night_pre = frac_night_pre
    this%frac_night_pet = frac_night_pet
    this%frac_night_temp = frac_night_temp
    this%frac_night_ssrd = frac_night_ssrd
    this%frac_night_strd = frac_night_strd
    this%share_frac = share_frac

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_meteo_from_file

  !> \brief Set config_meteo values
  integer function nml_config_meteo_set(this, &
    read_meteo_weights, &
    pre_weights_path, &
    pet_weights_path, &
    temp_weights_path, &
    ssrd_weights_path, &
    strd_weights_path, &
    pre_weights_var, &
    pet_weights_var, &
    temp_weights_var, &
    ssrd_weights_var, &
    strd_weights_var, &
    frac_night_pre, &
    frac_night_pet, &
    frac_night_temp, &
    frac_night_ssrd, &
    frac_night_strd, &
    share_frac, &
    errmsg) result(status)

    class(nml_config_meteo_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    logical, dimension(:), intent(in), optional :: read_meteo_weights !< Read meteorological weights
    character(len=*), dimension(:), intent(in), optional :: pre_weights_path !< Precipitation weights path
    character(len=*), dimension(:), intent(in), optional :: pet_weights_path !< Potential evapotranspiration weights path
    character(len=*), dimension(:), intent(in), optional :: temp_weights_path !< Surface downward shortwave radiation weights path
    character(len=*), dimension(:), intent(in), optional :: ssrd_weights_path !< Surface downward shortwave radiation weights path
    character(len=*), dimension(:), intent(in), optional :: strd_weights_path !< Surface downward longwave radiation weights path
    character(len=*), dimension(:), intent(in), optional :: pre_weights_var !< Precipitation weights variable name
    character(len=*), dimension(:), intent(in), optional :: pet_weights_var !< Potential evapotranspiration weights variable name
    character(len=*), dimension(:), intent(in), optional :: temp_weights_var !< Average temperature weights variable name
    character(len=*), dimension(:), intent(in), optional :: ssrd_weights_var !< Surface downward shortwave radiation weights variable name
    character(len=*), dimension(:), intent(in), optional :: strd_weights_var !< Surface downward longwave radiation weights variable name
    real(dp), dimension(:, :), intent(in), optional :: frac_night_pre !< Fraction of nightly precipitation
    real(dp), dimension(:, :), intent(in), optional :: frac_night_pet !< Fraction of nightly potential evapotranspiration
    real(dp), dimension(:, :), intent(in), optional :: frac_night_temp !< Fraction of nightly temperature
    real(dp), dimension(:, :), intent(in), optional :: frac_night_ssrd !< Fraction of nightly surface downward shortwave radiation
    real(dp), dimension(:, :), intent(in), optional :: frac_night_strd !< Fraction of nightly surface downward longwave radiation
    logical, intent(in), optional :: share_frac !< Share fractions between domains
    integer :: &
      lb__1, &
      lb__2, &
      ub__1, &
      ub__2

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(read_meteo_weights)) then
      if (size(read_meteo_weights, 1) > size(this%read_meteo_weights, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'read_meteo_weights'"
        return
      end if
      lb__1 = lbound(this%read_meteo_weights, 1)
      ub__1 = lb__1 + size(read_meteo_weights, 1) - 1
      this%read_meteo_weights(lb__1:ub__1) = read_meteo_weights
    end if
    if (present(pre_weights_path)) then
      if (size(pre_weights_path, 1) > size(this%pre_weights_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pre_weights_path'"
        return
      end if
      lb__1 = lbound(this%pre_weights_path, 1)
      ub__1 = lb__1 + size(pre_weights_path, 1) - 1
      this%pre_weights_path(lb__1:ub__1) = pre_weights_path
    end if
    if (present(pet_weights_path)) then
      if (size(pet_weights_path, 1) > size(this%pet_weights_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_weights_path'"
        return
      end if
      lb__1 = lbound(this%pet_weights_path, 1)
      ub__1 = lb__1 + size(pet_weights_path, 1) - 1
      this%pet_weights_path(lb__1:ub__1) = pet_weights_path
    end if
    if (present(temp_weights_path)) then
      if (size(temp_weights_path, 1) > size(this%temp_weights_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'temp_weights_path'"
        return
      end if
      lb__1 = lbound(this%temp_weights_path, 1)
      ub__1 = lb__1 + size(temp_weights_path, 1) - 1
      this%temp_weights_path(lb__1:ub__1) = temp_weights_path
    end if
    if (present(ssrd_weights_path)) then
      if (size(ssrd_weights_path, 1) > size(this%ssrd_weights_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ssrd_weights_path'"
        return
      end if
      lb__1 = lbound(this%ssrd_weights_path, 1)
      ub__1 = lb__1 + size(ssrd_weights_path, 1) - 1
      this%ssrd_weights_path(lb__1:ub__1) = ssrd_weights_path
    end if
    if (present(strd_weights_path)) then
      if (size(strd_weights_path, 1) > size(this%strd_weights_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'strd_weights_path'"
        return
      end if
      lb__1 = lbound(this%strd_weights_path, 1)
      ub__1 = lb__1 + size(strd_weights_path, 1) - 1
      this%strd_weights_path(lb__1:ub__1) = strd_weights_path
    end if
    if (present(pre_weights_var)) then
      if (size(pre_weights_var, 1) > size(this%pre_weights_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pre_weights_var'"
        return
      end if
      lb__1 = lbound(this%pre_weights_var, 1)
      ub__1 = lb__1 + size(pre_weights_var, 1) - 1
      this%pre_weights_var(lb__1:ub__1) = pre_weights_var
    end if
    if (present(pet_weights_var)) then
      if (size(pet_weights_var, 1) > size(this%pet_weights_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'pet_weights_var'"
        return
      end if
      lb__1 = lbound(this%pet_weights_var, 1)
      ub__1 = lb__1 + size(pet_weights_var, 1) - 1
      this%pet_weights_var(lb__1:ub__1) = pet_weights_var
    end if
    if (present(temp_weights_var)) then
      if (size(temp_weights_var, 1) > size(this%temp_weights_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'temp_weights_var'"
        return
      end if
      lb__1 = lbound(this%temp_weights_var, 1)
      ub__1 = lb__1 + size(temp_weights_var, 1) - 1
      this%temp_weights_var(lb__1:ub__1) = temp_weights_var
    end if
    if (present(ssrd_weights_var)) then
      if (size(ssrd_weights_var, 1) > size(this%ssrd_weights_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'ssrd_weights_var'"
        return
      end if
      lb__1 = lbound(this%ssrd_weights_var, 1)
      ub__1 = lb__1 + size(ssrd_weights_var, 1) - 1
      this%ssrd_weights_var(lb__1:ub__1) = ssrd_weights_var
    end if
    if (present(strd_weights_var)) then
      if (size(strd_weights_var, 1) > size(this%strd_weights_var, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'strd_weights_var'"
        return
      end if
      lb__1 = lbound(this%strd_weights_var, 1)
      ub__1 = lb__1 + size(strd_weights_var, 1) - 1
      this%strd_weights_var(lb__1:ub__1) = strd_weights_var
    end if
    if (present(frac_night_pre)) then
      if (size(frac_night_pre, 1) > size(this%frac_night_pre, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'frac_night_pre'"
        return
      end if
      lb__1 = lbound(this%frac_night_pre, 1)
      ub__1 = lb__1 + size(frac_night_pre, 1) - 1
      if (size(frac_night_pre, 2) > size(this%frac_night_pre, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'frac_night_pre'"
        return
      end if
      lb__2 = lbound(this%frac_night_pre, 2)
      ub__2 = lb__2 + size(frac_night_pre, 2) - 1
      this%frac_night_pre(lb__1:ub__1, lb__2:ub__2) = frac_night_pre
    end if
    if (present(frac_night_pet)) then
      if (size(frac_night_pet, 1) > size(this%frac_night_pet, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'frac_night_pet'"
        return
      end if
      lb__1 = lbound(this%frac_night_pet, 1)
      ub__1 = lb__1 + size(frac_night_pet, 1) - 1
      if (size(frac_night_pet, 2) > size(this%frac_night_pet, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'frac_night_pet'"
        return
      end if
      lb__2 = lbound(this%frac_night_pet, 2)
      ub__2 = lb__2 + size(frac_night_pet, 2) - 1
      this%frac_night_pet(lb__1:ub__1, lb__2:ub__2) = frac_night_pet
    end if
    if (present(frac_night_temp)) then
      if (size(frac_night_temp, 1) > size(this%frac_night_temp, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'frac_night_temp'"
        return
      end if
      lb__1 = lbound(this%frac_night_temp, 1)
      ub__1 = lb__1 + size(frac_night_temp, 1) - 1
      if (size(frac_night_temp, 2) > size(this%frac_night_temp, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'frac_night_temp'"
        return
      end if
      lb__2 = lbound(this%frac_night_temp, 2)
      ub__2 = lb__2 + size(frac_night_temp, 2) - 1
      this%frac_night_temp(lb__1:ub__1, lb__2:ub__2) = frac_night_temp
    end if
    if (present(frac_night_ssrd)) then
      if (size(frac_night_ssrd, 1) > size(this%frac_night_ssrd, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'frac_night_ssrd'"
        return
      end if
      lb__1 = lbound(this%frac_night_ssrd, 1)
      ub__1 = lb__1 + size(frac_night_ssrd, 1) - 1
      if (size(frac_night_ssrd, 2) > size(this%frac_night_ssrd, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'frac_night_ssrd'"
        return
      end if
      lb__2 = lbound(this%frac_night_ssrd, 2)
      ub__2 = lb__2 + size(frac_night_ssrd, 2) - 1
      this%frac_night_ssrd(lb__1:ub__1, lb__2:ub__2) = frac_night_ssrd
    end if
    if (present(frac_night_strd)) then
      if (size(frac_night_strd, 1) > size(this%frac_night_strd, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'frac_night_strd'"
        return
      end if
      lb__1 = lbound(this%frac_night_strd, 1)
      ub__1 = lb__1 + size(frac_night_strd, 1) - 1
      if (size(frac_night_strd, 2) > size(this%frac_night_strd, 2)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 2 exceeds bounds for 'frac_night_strd'"
        return
      end if
      lb__2 = lbound(this%frac_night_strd, 2)
      ub__2 = lb__2 + size(frac_night_strd, 2) - 1
      this%frac_night_strd(lb__1:ub__1, lb__2:ub__2) = frac_night_strd
    end if
    if (present(share_frac)) this%share_frac = share_frac

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_meteo_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_meteo_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_meteo_t), intent(in) :: this !< namelist instance
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
    case ("read_meteo_weights")
      if (.not. allocated(this%read_meteo_weights)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%read_meteo_weights), ubound(this%read_meteo_weights), &
          "read_meteo_weights", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("pre_weights_path")
      if (.not. allocated(this%pre_weights_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pre_weights_path), ubound(this%pre_weights_path), &
          "pre_weights_path", errmsg)
        if (status /= NML_OK) return
        if (this%pre_weights_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%pre_weights_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("pet_weights_path")
      if (.not. allocated(this%pet_weights_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_weights_path), ubound(this%pet_weights_path), &
          "pet_weights_path", errmsg)
        if (status /= NML_OK) return
        if (this%pet_weights_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%pet_weights_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("temp_weights_path")
      if (.not. allocated(this%temp_weights_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%temp_weights_path), ubound(this%temp_weights_path), &
          "temp_weights_path", errmsg)
        if (status /= NML_OK) return
        if (this%temp_weights_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%temp_weights_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("ssrd_weights_path")
      if (.not. allocated(this%ssrd_weights_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ssrd_weights_path), ubound(this%ssrd_weights_path), &
          "ssrd_weights_path", errmsg)
        if (status /= NML_OK) return
        if (this%ssrd_weights_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%ssrd_weights_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("strd_weights_path")
      if (.not. allocated(this%strd_weights_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%strd_weights_path), ubound(this%strd_weights_path), &
          "strd_weights_path", errmsg)
        if (status /= NML_OK) return
        if (this%strd_weights_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%strd_weights_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("pre_weights_var")
      if (.not. allocated(this%pre_weights_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pre_weights_var), ubound(this%pre_weights_var), &
          "pre_weights_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("pet_weights_var")
      if (.not. allocated(this%pet_weights_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%pet_weights_var), ubound(this%pet_weights_var), &
          "pet_weights_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("temp_weights_var")
      if (.not. allocated(this%temp_weights_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%temp_weights_var), ubound(this%temp_weights_var), &
          "temp_weights_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("ssrd_weights_var")
      if (.not. allocated(this%ssrd_weights_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%ssrd_weights_var), ubound(this%ssrd_weights_var), &
          "ssrd_weights_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("strd_weights_var")
      if (.not. allocated(this%strd_weights_var)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%strd_weights_var), ubound(this%strd_weights_var), &
          "strd_weights_var", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("frac_night_pre")
      if (.not. allocated(this%frac_night_pre)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%frac_night_pre), ubound(this%frac_night_pre), &
          "frac_night_pre", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%frac_night_pre(idx(1), idx(2)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%frac_night_pre))) status = NML_ERR_NOT_SET
      end if
    case ("frac_night_pet")
      if (.not. allocated(this%frac_night_pet)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%frac_night_pet), ubound(this%frac_night_pet), &
          "frac_night_pet", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%frac_night_pet(idx(1), idx(2)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%frac_night_pet))) status = NML_ERR_NOT_SET
      end if
    case ("frac_night_temp")
      if (.not. allocated(this%frac_night_temp)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%frac_night_temp), ubound(this%frac_night_temp), &
          "frac_night_temp", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%frac_night_temp(idx(1), idx(2)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%frac_night_temp))) status = NML_ERR_NOT_SET
      end if
    case ("frac_night_ssrd")
      if (.not. allocated(this%frac_night_ssrd)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%frac_night_ssrd), ubound(this%frac_night_ssrd), &
          "frac_night_ssrd", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%frac_night_ssrd(idx(1), idx(2)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%frac_night_ssrd))) status = NML_ERR_NOT_SET
      end if
    case ("frac_night_strd")
      if (.not. allocated(this%frac_night_strd)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%frac_night_strd), ubound(this%frac_night_strd), &
          "frac_night_strd", errmsg)
        if (status /= NML_OK) return
        if (ieee_is_nan(this%frac_night_strd(idx(1), idx(2)))) status = NML_ERR_NOT_SET
      else
        if (all(ieee_is_nan(this%frac_night_strd))) status = NML_ERR_NOT_SET
      end if
    case ("share_frac")
      if (present(idx)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "index not supported for 'share_frac'"
        return
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_meteo_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_meteo_is_valid(this, errmsg) result(status)
    class(nml_config_meteo_t), intent(in) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    integer :: istat

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    if (.not. this%is_configured) then
      status = NML_ERR_NOT_SET
      if (present(errmsg)) errmsg = "namelist not configured; call set or from_file"
      return
    end if

  end function nml_config_meteo_is_valid

end module nml_config_meteo
