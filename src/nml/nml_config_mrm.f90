!> \file nml_config_mrm.f90
!> \copydoc nml_config_mrm

!> \brief mRM configuration
!> \details Configuration for the multi-scale routing model (mRM) in mHM.
!> \version 0.2
!> \authors Sebastian Mueller
!> \date    Jun 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_namelists
module nml_config_mrm
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
    i4, &
    dp

  implicit none

  ! default values
  logical, parameter, public :: river_net_order_root_based__default = .false.
  integer(i4), parameter, public :: river_net_omp_level_min__default = -1_i4
  integer(i4), parameter, public :: max_route_step__default = 86400_i4
  integer(i4), parameter, public :: upscale_mode__default = 1_i4
  real(dp), parameter, public :: length_percentile__default = 40.0_dp
  logical, parameter, public :: read_restart__default = .false.
  logical, parameter, public :: read_restart_fluxes__default = .true.
  logical, parameter, public :: write_restart__default = .false.

  ! enum values
  integer(i4), parameter, public :: max_route_step__enum_values(19) = [60_i4, 120_i4, 180_i4, 240_i4, 300_i4, 360_i4, 600_i4, 720_i4, 900_i4, 1200_i4, 1800_i4, 3600_i4, 7200_i4, 10800_i4, 14400_i4, 21600_i4, 28800_i4, 43200_i4, 86400_i4]
  integer(i4), parameter, public :: upscale_mode__enum_values(2) = [0_i4, 1_i4]

  ! bounds values
  integer(i4), parameter, public :: river_net_omp_level_min__min = -1_i4
  real(dp), parameter, public :: length_percentile__min = 0.0_dp
  real(dp), parameter, public :: length_percentile__max = 100.0_dp

  !> \class nml_config_mrm_t
  !> \brief mRM configuration
  !> \details Configuration for the multi-scale routing model (mRM) in mHM.
  type, public :: nml_config_mrm_t
    logical :: is_configured = .false. !< whether the namelist has been configured
    integer :: n_domains = n_domains__default !< runtime dimension for n_domains
    logical, allocatable, dimension(:) :: river_net_order_root_based !< Flag for root based river network ordering.
    integer(i4), allocatable, dimension(:) :: river_net_omp_level_min !< Minimum level size for OpenMP parallelization.
    integer(i4), allocatable, dimension(:) :: max_route_step !< Maximum numerical routing substep in seconds.
    integer(i4), allocatable, dimension(:) :: upscale_mode !< River upscaling mode.
    real(dp), allocatable, dimension(:) :: length_percentile !< Percentile for the minimum upscaled link length.
    character(len=buf), allocatable, dimension(:) :: scc_gauges_path !< Path for SCC gauges NetCDF file.
    character(len=buf), allocatable, dimension(:) :: output_path !< Path for output file.
    character(len=buf), allocatable, dimension(:) :: output_node_path !< Path for node based output file.
    logical, allocatable, dimension(:) :: read_restart !< Read restart
    logical, allocatable, dimension(:) :: read_restart_fluxes !< Read restart fluxes
    character(len=buf), allocatable, dimension(:) :: restart_input_path !< Restart input path
    logical, allocatable, dimension(:) :: write_restart !< Write restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path !< Restart output path
    character(len=buf), allocatable, dimension(:) :: diagnostics_path !< Diagnostics output path
  contains
    procedure :: init => nml_config_mrm_init
    procedure :: set_dims => nml_config_mrm_set_dims
    procedure :: from_file => nml_config_mrm_from_file
    procedure :: set => nml_config_mrm_set
    procedure :: is_set => nml_config_mrm_is_set
    procedure :: is_valid => nml_config_mrm_is_valid
  end type nml_config_mrm_t

contains

  !> \brief Check whether a value is part of an enum
  elemental logical function max_route_step__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == max_route_step__enum_values)
  end function max_route_step__in_enum

  !> \brief Check whether a value is part of an enum
  elemental logical function upscale_mode__in_enum(val, allow_missing) result(in_enum)
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
    in_enum = any(val == upscale_mode__enum_values)
  end function upscale_mode__in_enum

  !> \brief Check whether a value is within bounds
  elemental logical function river_net_omp_level_min__in_bounds(val, allow_missing) result(in_bounds)
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
    if (val < river_net_omp_level_min__min) in_bounds = .false.
  end function river_net_omp_level_min__in_bounds

  !> \brief Check whether a value is within bounds
  elemental logical function length_percentile__in_bounds(val, allow_missing) result(in_bounds)
    real(dp), intent(in) :: val !< value to check
    logical, intent(in), optional :: allow_missing !< allow sentinel values as valid

    if (present(allow_missing)) then
      if (allow_missing) then
        if (ieee_is_nan(val)) then
          in_bounds = .true.
          return
        end if
      end if
    end if

    in_bounds = .true.
    if (val < length_percentile__min) in_bounds = .false.
    if (val > length_percentile__max) in_bounds = .false.
  end function length_percentile__in_bounds

  !> \brief Initialize defaults and sentinels for config_mrm
  integer function nml_config_mrm_init(this, errmsg) result(status)
    class(nml_config_mrm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values

    status = NML_OK
    if (present(errmsg)) errmsg = ""
    this%is_configured = .false.

    ! allocate runtime-sized fields
    if (allocated(this%river_net_order_root_based)) deallocate(this%river_net_order_root_based)
    allocate(this%river_net_order_root_based(this%n_domains))
    if (allocated(this%river_net_omp_level_min)) deallocate(this%river_net_omp_level_min)
    allocate(this%river_net_omp_level_min(this%n_domains))
    if (allocated(this%max_route_step)) deallocate(this%max_route_step)
    allocate(this%max_route_step(this%n_domains))
    if (allocated(this%upscale_mode)) deallocate(this%upscale_mode)
    allocate(this%upscale_mode(this%n_domains))
    if (allocated(this%length_percentile)) deallocate(this%length_percentile)
    allocate(this%length_percentile(this%n_domains))
    if (allocated(this%scc_gauges_path)) deallocate(this%scc_gauges_path)
    allocate(character(len=buf) :: this%scc_gauges_path(this%n_domains))
    if (allocated(this%output_path)) deallocate(this%output_path)
    allocate(character(len=buf) :: this%output_path(this%n_domains))
    if (allocated(this%output_node_path)) deallocate(this%output_node_path)
    allocate(character(len=buf) :: this%output_node_path(this%n_domains))
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    allocate(this%read_restart(this%n_domains))
    if (allocated(this%read_restart_fluxes)) deallocate(this%read_restart_fluxes)
    allocate(this%read_restart_fluxes(this%n_domains))
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    allocate(character(len=buf) :: this%restart_input_path(this%n_domains))
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    allocate(this%write_restart(this%n_domains))
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    allocate(character(len=buf) :: this%restart_output_path(this%n_domains))
    if (allocated(this%diagnostics_path)) deallocate(this%diagnostics_path)
    allocate(character(len=buf) :: this%diagnostics_path(this%n_domains))

    ! sentinel values for required/optional parameters
    this%scc_gauges_path = achar(0) ! sentinel for optional string array
    this%output_path = achar(0) ! sentinel for optional string array
    this%output_node_path = achar(0) ! sentinel for optional string array
    this%restart_input_path = achar(0) ! sentinel for optional string array
    this%restart_output_path = achar(0) ! sentinel for optional string array
    this%diagnostics_path = achar(0) ! sentinel for optional string array
    ! default values
    this%river_net_order_root_based = river_net_order_root_based__default
    this%river_net_omp_level_min = river_net_omp_level_min__default
    this%max_route_step = max_route_step__default
    this%upscale_mode = upscale_mode__default
    this%length_percentile = length_percentile__default
    this%read_restart = read_restart__default
    this%read_restart_fluxes = read_restart_fluxes__default
    this%write_restart = write_restart__default
  end function nml_config_mrm_init

  !> \brief Reset runtime dimensions for config_mrm
  integer function nml_config_mrm_set_dims(this, &
    n_domains, &
    errmsg) result(status)
    class(nml_config_mrm_t), intent(inout) :: this !< namelist instance
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
    if (allocated(this%river_net_order_root_based)) deallocate(this%river_net_order_root_based)
    if (allocated(this%river_net_omp_level_min)) deallocate(this%river_net_omp_level_min)
    if (allocated(this%max_route_step)) deallocate(this%max_route_step)
    if (allocated(this%upscale_mode)) deallocate(this%upscale_mode)
    if (allocated(this%length_percentile)) deallocate(this%length_percentile)
    if (allocated(this%scc_gauges_path)) deallocate(this%scc_gauges_path)
    if (allocated(this%output_path)) deallocate(this%output_path)
    if (allocated(this%output_node_path)) deallocate(this%output_node_path)
    if (allocated(this%read_restart)) deallocate(this%read_restart)
    if (allocated(this%read_restart_fluxes)) deallocate(this%read_restart_fluxes)
    if (allocated(this%restart_input_path)) deallocate(this%restart_input_path)
    if (allocated(this%write_restart)) deallocate(this%write_restart)
    if (allocated(this%restart_output_path)) deallocate(this%restart_output_path)
    if (allocated(this%diagnostics_path)) deallocate(this%diagnostics_path)
    this%is_configured = .false.
  end function nml_config_mrm_set_dims


  !> \brief Read config_mrm namelist from file
  integer function nml_config_mrm_from_file(this, file, errmsg) result(status)
    class(nml_config_mrm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(in) :: file !< path to namelist file
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    ! namelist variables
    logical, allocatable, dimension(:) :: river_net_order_root_based
    integer(i4), allocatable, dimension(:) :: river_net_omp_level_min
    integer(i4), allocatable, dimension(:) :: max_route_step
    integer(i4), allocatable, dimension(:) :: upscale_mode
    real(dp), allocatable, dimension(:) :: length_percentile
    character(len=buf), allocatable, dimension(:) :: scc_gauges_path
    character(len=buf), allocatable, dimension(:) :: output_path
    character(len=buf), allocatable, dimension(:) :: output_node_path
    logical, allocatable, dimension(:) :: read_restart
    logical, allocatable, dimension(:) :: read_restart_fluxes
    character(len=buf), allocatable, dimension(:) :: restart_input_path
    logical, allocatable, dimension(:) :: write_restart
    character(len=buf), allocatable, dimension(:) :: restart_output_path
    character(len=buf), allocatable, dimension(:) :: diagnostics_path
    ! locals
    type(nml_file_t) :: nml
    integer :: iostat
    integer :: close_status
    character(len=nml_line_buffer) :: iomsg

    namelist /config_mrm/ &
      river_net_order_root_based, &
      river_net_omp_level_min, &
      max_route_step, &
      upscale_mode, &
      length_percentile, &
      scc_gauges_path, &
      output_path, &
      output_node_path, &
      read_restart, &
      read_restart_fluxes, &
      restart_input_path, &
      write_restart, &
      restart_output_path, &
      diagnostics_path

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return
    ! allocate local namelist variables matching runtime-sized fields
    if (allocated(river_net_order_root_based)) deallocate(river_net_order_root_based)
    allocate(river_net_order_root_based(this%n_domains))
    if (allocated(river_net_omp_level_min)) deallocate(river_net_omp_level_min)
    allocate(river_net_omp_level_min(this%n_domains))
    if (allocated(max_route_step)) deallocate(max_route_step)
    allocate(max_route_step(this%n_domains))
    if (allocated(upscale_mode)) deallocate(upscale_mode)
    allocate(upscale_mode(this%n_domains))
    if (allocated(length_percentile)) deallocate(length_percentile)
    allocate(length_percentile(this%n_domains))
    if (allocated(scc_gauges_path)) deallocate(scc_gauges_path)
    allocate(character(len=buf) :: scc_gauges_path(this%n_domains))
    if (allocated(output_path)) deallocate(output_path)
    allocate(character(len=buf) :: output_path(this%n_domains))
    if (allocated(output_node_path)) deallocate(output_node_path)
    allocate(character(len=buf) :: output_node_path(this%n_domains))
    if (allocated(read_restart)) deallocate(read_restart)
    allocate(read_restart(this%n_domains))
    if (allocated(read_restart_fluxes)) deallocate(read_restart_fluxes)
    allocate(read_restart_fluxes(this%n_domains))
    if (allocated(restart_input_path)) deallocate(restart_input_path)
    allocate(character(len=buf) :: restart_input_path(this%n_domains))
    if (allocated(write_restart)) deallocate(write_restart)
    allocate(write_restart(this%n_domains))
    if (allocated(restart_output_path)) deallocate(restart_output_path)
    allocate(character(len=buf) :: restart_output_path(this%n_domains))
    if (allocated(diagnostics_path)) deallocate(diagnostics_path)
    allocate(character(len=buf) :: diagnostics_path(this%n_domains))
    river_net_order_root_based = this%river_net_order_root_based
    river_net_omp_level_min = this%river_net_omp_level_min
    max_route_step = this%max_route_step
    upscale_mode = this%upscale_mode
    length_percentile = this%length_percentile
    scc_gauges_path = this%scc_gauges_path
    output_path = this%output_path
    output_node_path = this%output_node_path
    read_restart = this%read_restart
    read_restart_fluxes = this%read_restart_fluxes
    restart_input_path = this%restart_input_path
    write_restart = this%write_restart
    restart_output_path = this%restart_output_path
    diagnostics_path = this%diagnostics_path

    status = nml%open(file, errmsg=errmsg)
    if (status /= NML_OK) return

    status = nml%find("config_mrm", errmsg=errmsg)
    if (status /= NML_OK) then
      close_status = nml%close()
      return
    end if

    ! read namelist
    read(nml%unit, nml=config_mrm, iostat=iostat, iomsg=iomsg)
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
    this%river_net_order_root_based = river_net_order_root_based
    this%river_net_omp_level_min = river_net_omp_level_min
    this%max_route_step = max_route_step
    this%upscale_mode = upscale_mode
    this%length_percentile = length_percentile
    this%scc_gauges_path = scc_gauges_path
    this%output_path = output_path
    this%output_node_path = output_node_path
    this%read_restart = read_restart
    this%read_restart_fluxes = read_restart_fluxes
    this%restart_input_path = restart_input_path
    this%write_restart = write_restart
    this%restart_output_path = restart_output_path
    this%diagnostics_path = diagnostics_path

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mrm_from_file

  !> \brief Set config_mrm values
  integer function nml_config_mrm_set(this, &
    river_net_order_root_based, &
    river_net_omp_level_min, &
    max_route_step, &
    upscale_mode, &
    length_percentile, &
    scc_gauges_path, &
    output_path, &
    output_node_path, &
    read_restart, &
    read_restart_fluxes, &
    restart_input_path, &
    write_restart, &
    restart_output_path, &
    diagnostics_path, &
    errmsg) result(status)

    class(nml_config_mrm_t), intent(inout) :: this !< namelist instance
    character(len=*), intent(out), optional :: errmsg !< error message for non-OK status values
    logical, dimension(:), intent(in), optional :: river_net_order_root_based !< Flag for root based river network ordering.
    integer(i4), dimension(:), intent(in), optional :: river_net_omp_level_min !< Minimum level size for OpenMP parallelization.
    integer(i4), dimension(:), intent(in), optional :: max_route_step !< Maximum numerical routing substep in seconds.
    integer(i4), dimension(:), intent(in), optional :: upscale_mode !< River upscaling mode.
    real(dp), dimension(:), intent(in), optional :: length_percentile !< Percentile for the minimum upscaled link length.
    character(len=*), dimension(:), intent(in), optional :: scc_gauges_path !< Path for SCC gauges NetCDF file.
    character(len=*), dimension(:), intent(in), optional :: output_path !< Path for output file.
    character(len=*), dimension(:), intent(in), optional :: output_node_path !< Path for node based output file.
    logical, dimension(:), intent(in), optional :: read_restart !< Read restart
    logical, dimension(:), intent(in), optional :: read_restart_fluxes !< Read restart fluxes
    character(len=*), dimension(:), intent(in), optional :: restart_input_path !< Restart input path
    logical, dimension(:), intent(in), optional :: write_restart !< Write restart
    character(len=*), dimension(:), intent(in), optional :: restart_output_path !< Restart output path
    character(len=*), dimension(:), intent(in), optional :: diagnostics_path !< Diagnostics output path
    integer :: &
      lb__1, &
      ub__1

    status = this%init(errmsg=errmsg)
    if (status /= NML_OK) return

    ! required parameters
    ! override with provided values
    if (present(river_net_order_root_based)) then
      if (size(river_net_order_root_based, 1) > size(this%river_net_order_root_based, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'river_net_order_root_based'"
        return
      end if
      lb__1 = lbound(this%river_net_order_root_based, 1)
      ub__1 = lb__1 + size(river_net_order_root_based, 1) - 1
      this%river_net_order_root_based(lb__1:ub__1) = river_net_order_root_based
    end if
    if (present(river_net_omp_level_min)) then
      if (size(river_net_omp_level_min, 1) > size(this%river_net_omp_level_min, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'river_net_omp_level_min'"
        return
      end if
      lb__1 = lbound(this%river_net_omp_level_min, 1)
      ub__1 = lb__1 + size(river_net_omp_level_min, 1) - 1
      this%river_net_omp_level_min(lb__1:ub__1) = river_net_omp_level_min
    end if
    if (present(max_route_step)) then
      if (size(max_route_step, 1) > size(this%max_route_step, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'max_route_step'"
        return
      end if
      lb__1 = lbound(this%max_route_step, 1)
      ub__1 = lb__1 + size(max_route_step, 1) - 1
      this%max_route_step(lb__1:ub__1) = max_route_step
    end if
    if (present(upscale_mode)) then
      if (size(upscale_mode, 1) > size(this%upscale_mode, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'upscale_mode'"
        return
      end if
      lb__1 = lbound(this%upscale_mode, 1)
      ub__1 = lb__1 + size(upscale_mode, 1) - 1
      this%upscale_mode(lb__1:ub__1) = upscale_mode
    end if
    if (present(length_percentile)) then
      if (size(length_percentile, 1) > size(this%length_percentile, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'length_percentile'"
        return
      end if
      lb__1 = lbound(this%length_percentile, 1)
      ub__1 = lb__1 + size(length_percentile, 1) - 1
      this%length_percentile(lb__1:ub__1) = length_percentile
    end if
    if (present(scc_gauges_path)) then
      if (size(scc_gauges_path, 1) > size(this%scc_gauges_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'scc_gauges_path'"
        return
      end if
      lb__1 = lbound(this%scc_gauges_path, 1)
      ub__1 = lb__1 + size(scc_gauges_path, 1) - 1
      this%scc_gauges_path(lb__1:ub__1) = scc_gauges_path
    end if
    if (present(output_path)) then
      if (size(output_path, 1) > size(this%output_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'output_path'"
        return
      end if
      lb__1 = lbound(this%output_path, 1)
      ub__1 = lb__1 + size(output_path, 1) - 1
      this%output_path(lb__1:ub__1) = output_path
    end if
    if (present(output_node_path)) then
      if (size(output_node_path, 1) > size(this%output_node_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'output_node_path'"
        return
      end if
      lb__1 = lbound(this%output_node_path, 1)
      ub__1 = lb__1 + size(output_node_path, 1) - 1
      this%output_node_path(lb__1:ub__1) = output_node_path
    end if
    if (present(read_restart)) then
      if (size(read_restart, 1) > size(this%read_restart, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'read_restart'"
        return
      end if
      lb__1 = lbound(this%read_restart, 1)
      ub__1 = lb__1 + size(read_restart, 1) - 1
      this%read_restart(lb__1:ub__1) = read_restart
    end if
    if (present(read_restart_fluxes)) then
      if (size(read_restart_fluxes, 1) > size(this%read_restart_fluxes, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'read_restart_fluxes'"
        return
      end if
      lb__1 = lbound(this%read_restart_fluxes, 1)
      ub__1 = lb__1 + size(read_restart_fluxes, 1) - 1
      this%read_restart_fluxes(lb__1:ub__1) = read_restart_fluxes
    end if
    if (present(restart_input_path)) then
      if (size(restart_input_path, 1) > size(this%restart_input_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'restart_input_path'"
        return
      end if
      lb__1 = lbound(this%restart_input_path, 1)
      ub__1 = lb__1 + size(restart_input_path, 1) - 1
      this%restart_input_path(lb__1:ub__1) = restart_input_path
    end if
    if (present(write_restart)) then
      if (size(write_restart, 1) > size(this%write_restart, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'write_restart'"
        return
      end if
      lb__1 = lbound(this%write_restart, 1)
      ub__1 = lb__1 + size(write_restart, 1) - 1
      this%write_restart(lb__1:ub__1) = write_restart
    end if
    if (present(restart_output_path)) then
      if (size(restart_output_path, 1) > size(this%restart_output_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'restart_output_path'"
        return
      end if
      lb__1 = lbound(this%restart_output_path, 1)
      ub__1 = lb__1 + size(restart_output_path, 1) - 1
      this%restart_output_path(lb__1:ub__1) = restart_output_path
    end if
    if (present(diagnostics_path)) then
      if (size(diagnostics_path, 1) > size(this%diagnostics_path, 1)) then
        status = NML_ERR_INVALID_INDEX
        if (present(errmsg)) errmsg = "dimension 1 exceeds bounds for 'diagnostics_path'"
        return
      end if
      lb__1 = lbound(this%diagnostics_path, 1)
      ub__1 = lb__1 + size(diagnostics_path, 1) - 1
      this%diagnostics_path(lb__1:ub__1) = diagnostics_path
    end if

    ! mark as configured
    this%is_configured = .true.
    status = NML_OK
  end function nml_config_mrm_set

  !> \brief Check whether a namelist value was set
  integer function nml_config_mrm_is_set(this, name, idx, errmsg) result(status)
    class(nml_config_mrm_t), intent(in) :: this !< namelist instance
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
    case ("river_net_order_root_based")
      if (.not. allocated(this%river_net_order_root_based)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%river_net_order_root_based), ubound(this%river_net_order_root_based), &
          "river_net_order_root_based", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("river_net_omp_level_min")
      if (.not. allocated(this%river_net_omp_level_min)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%river_net_omp_level_min), ubound(this%river_net_omp_level_min), &
          "river_net_omp_level_min", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("max_route_step")
      if (.not. allocated(this%max_route_step)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%max_route_step), ubound(this%max_route_step), &
          "max_route_step", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("upscale_mode")
      if (.not. allocated(this%upscale_mode)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%upscale_mode), ubound(this%upscale_mode), &
          "upscale_mode", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("length_percentile")
      if (.not. allocated(this%length_percentile)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%length_percentile), ubound(this%length_percentile), &
          "length_percentile", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("scc_gauges_path")
      if (.not. allocated(this%scc_gauges_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%scc_gauges_path), ubound(this%scc_gauges_path), &
          "scc_gauges_path", errmsg)
        if (status /= NML_OK) return
        if (this%scc_gauges_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%scc_gauges_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("output_path")
      if (.not. allocated(this%output_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%output_path), ubound(this%output_path), &
          "output_path", errmsg)
        if (status /= NML_OK) return
        if (this%output_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%output_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("output_node_path")
      if (.not. allocated(this%output_node_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%output_node_path), ubound(this%output_node_path), &
          "output_node_path", errmsg)
        if (status /= NML_OK) return
        if (this%output_node_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%output_node_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("read_restart")
      if (.not. allocated(this%read_restart)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%read_restart), ubound(this%read_restart), &
          "read_restart", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("read_restart_fluxes")
      if (.not. allocated(this%read_restart_fluxes)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%read_restart_fluxes), ubound(this%read_restart_fluxes), &
          "read_restart_fluxes", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("restart_input_path")
      if (.not. allocated(this%restart_input_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%restart_input_path), ubound(this%restart_input_path), &
          "restart_input_path", errmsg)
        if (status /= NML_OK) return
        if (this%restart_input_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%restart_input_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("write_restart")
      if (.not. allocated(this%write_restart)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%write_restart), ubound(this%write_restart), &
          "write_restart", errmsg)
        if (status /= NML_OK) return
      else
      end if
    case ("restart_output_path")
      if (.not. allocated(this%restart_output_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%restart_output_path), ubound(this%restart_output_path), &
          "restart_output_path", errmsg)
        if (status /= NML_OK) return
        if (this%restart_output_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%restart_output_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case ("diagnostics_path")
      if (.not. allocated(this%diagnostics_path)) then
        status = NML_ERR_NOT_SET
        return
      end if
      if (present(idx)) then
        status = idx_check(idx, lbound(this%diagnostics_path), ubound(this%diagnostics_path), &
          "diagnostics_path", errmsg)
        if (status /= NML_OK) return
        if (this%diagnostics_path(idx(1)) == achar(0)) status = NML_ERR_NOT_SET
      else
        if (all(this%diagnostics_path == achar(0))) status = NML_ERR_NOT_SET
      end if
    case default
      status = NML_ERR_INVALID_NAME
      if (present(errmsg)) errmsg = "unknown field: " // trim(name)
    end select
    if (status == NML_ERR_NOT_SET .and. present(errmsg)) then
      if (len_trim(errmsg) == 0) errmsg = "field not set: " // trim(name)
    end if
  end function nml_config_mrm_is_set

  !> \brief Validate required values and constraints
  integer function nml_config_mrm_is_valid(this, errmsg) result(status)
    class(nml_config_mrm_t), intent(in) :: this !< namelist instance
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
    if (allocated(this%max_route_step)) then
    if (.not. all(max_route_step__in_enum(this%max_route_step, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: max_route_step"
      return
    end if
    end if
    if (allocated(this%upscale_mode)) then
    if (.not. all(upscale_mode__in_enum(this%upscale_mode, allow_missing=.true.))) then
      status = NML_ERR_ENUM
      if (present(errmsg)) errmsg = "enum constraint failed: upscale_mode"
      return
    end if
    end if
    ! bounds constraints
    if (allocated(this%river_net_omp_level_min)) then
    if (.not. all(river_net_omp_level_min__in_bounds(this%river_net_omp_level_min, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: river_net_omp_level_min"
      return
    end if
    end if
    if (allocated(this%length_percentile)) then
    if (.not. all(length_percentile__in_bounds(this%length_percentile, allow_missing=.true.))) then
      status = NML_ERR_BOUNDS
      if (present(errmsg)) errmsg = "bounds constraint failed: length_percentile"
      return
    end if
    end if
  end function nml_config_mrm_is_valid

end module nml_config_mrm
