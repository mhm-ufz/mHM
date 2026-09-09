!> \file mo_mlm_container.f90
!> \brief Minimal v6 lake component.
!> \details Provides the process--1 delayed pass-through contract between mRM updates.
!> \changelog
!! - Pallav Shrestha (2018-2023): v5 lake and SCC coupling concepts.
!! - Sebastian Mueller (2026): v6 domain-owned mLM pass-through component.
!> \authors Pallav Shrestha
!> \authors Sebastian Mueller
!> \ingroup f_exchange
#include "logging.h"
module mo_mlm_container
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_logging
  use mo_kind, only: i4, i8, dp
  use mo_datetime, only: datetime, one_hour
  use mo_exchange_type, only: exchange_t
  use mo_grid, only: grid_t
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable
  use mo_lake_balance, only: lake_balance_outflow
  use mo_lake_topology, only: lake_derive_area, lake_read_topology, lake_write_topology
  use mo_grid_io, only: var, daily, monthly, yearly, no_time, time_units_delta
  use mo_points, only: points_t, cartesian, spherical
  use mo_points_io, only: points_output_dataset
  use mo_utils, only: is_close
  use mo_string_utils, only: n2s => num2str
  use nml_config_mlm, only: nml_config_mlm_t
  use nml_output_mlm, only: nml_output_mlm_t
  use nml_helper, only: NML_OK
  implicit none

  character(*), parameter :: s = "mlm" !< logging scope

  !> \class mlm_t
  !> \brief Domain-owned lake pass-through component.
  !> \authors Pallav Shrestha
  !> \authors Sebastian Mueller
  type, public :: mlm_t
    type(nml_config_mlm_t) :: config !< mLM configuration
    type(nml_output_mlm_t) :: output_config !< lake-point output configuration
    type(exchange_t), pointer :: exchange => null() !< owning domain exchange
    integer(i8), allocatable :: lake_ids(:) !< stable IDs in point-set order
    integer(i8), allocatable :: lake_map(:) !< stable IDs on packed level-0 lake cells
    real(dp), allocatable :: lake_area(:) !< lake surface areas [m2]
    type(grid_t) :: static_lake_grid !< lake footprint geometry retained for restart
    type(points_t) :: static_lake_points !< static point metadata retained for restart output
    real(dp), allocatable :: static_lake_max_levels(:) !< maximum levels retained for restart output
    type(points_t) :: restart_lake_points !< point set restored from restart metadata
    integer(i8), allocatable :: restart_lake_ids(:) !< stable IDs in restart-file order
    real(dp), allocatable :: restart_lake_max_levels(:) !< maximum levels restored from restart metadata
    real(dp), allocatable :: outflow(:) !< current lake outflow [m3 s-1]
    logical :: active = .false. !< whether process -1 is selected
    logical :: read_restart = .false. !< read mLM restart
    logical :: write_restart = .false. !< write mLM restart
    logical :: output_active = .false. !< whether lake-point output is enabled
    logical :: owns_restart_lake_points = .false. !< whether mLM published the restart point set
    logical :: owns_restart_lake_metadata = .false. !< whether mLM published restart IDs and maximum levels
    logical :: owns_restart_lake_grid = .false. !< whether mLM restored and published lake topology
    character(:), allocatable :: restart_input_path !< resolved input restart path
    character(:), allocatable :: restart_output_path !< resolved output restart path
    character(:), allocatable :: output_path !< resolved lake-point output path
    type(points_output_dataset) :: ds_out !< lake-point output dataset
  contains
    procedure :: set_dims => mlm_set_dims
    procedure :: configure => mlm_configure
    procedure :: connect => mlm_connect
    procedure :: initialize => mlm_initialize
    procedure :: update => mlm_update
    procedure :: finalize => mlm_finalize
    procedure, private :: create_restart => mlm_create_restart
    procedure, private :: read_restart_metadata => mlm_read_restart_metadata
    procedure, private :: read_restart_topology => mlm_read_restart_topology
    procedure, private :: read_restart_state => mlm_read_restart_state
    procedure, private :: validate_restart_metadata => mlm_validate_restart_metadata
    procedure, private :: create_output => mlm_create_output
    procedure, private :: validate_output_timing => mlm_validate_output_timing
    procedure, private :: at_output_boundary => mlm_at_output_boundary
    procedure, private :: update_output => mlm_update_output
  end type mlm_t

contains

  !> \brief Set generated lake-configuration dimensions.
  subroutine mlm_set_dims(self)
    class(mlm_t), target, intent(inout) :: self
    character(1024) :: errmsg
    integer :: status
    status = self%config%set_dims(n_domains=self%exchange%nml_n_domains, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "mLM: error setting config_mlm dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine mlm_set_dims

  !> \brief Configure process -1, restart paths, and optional lake-point output.
  subroutine mlm_configure(self, file, out_file)
    class(mlm_t), target, intent(inout) :: self
    character(*), optional, intent(in) :: file
    character(*), optional, intent(in) :: out_file
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer(i4) :: id(1), lake_case
    integer :: status

    lake_case = self%exchange%config%processes%lake
    self%active = lake_case /= 0_i4
    self%exchange%lake_pre%required = self%exchange%lake_pre%required .or. lake_case == -2_i4
    self%exchange%lake_pet%required = self%exchange%lake_pet%required .or. lake_case == -2_i4
    if (lake_case /= 0_i4 .and. lake_case /= -1_i4 .and. lake_case /= -2_i4) then
      log_fatal(*) "mLM: unsupported lake process case: ", n2s(lake_case)
      error stop 1
    end if
    if (.not.self%active) return
    if (self%exchange%config%processes%routing == 0_i4) then
      log_fatal(*) "mLM lake processes require active mRM routing."
      error stop 1
    end if
    id(1) = self%exchange%nml_domain_id
    if (present(file)) then
      path = self%exchange%get_path(file)
      status = self%config%from_file(path, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mLM: error reading config_mlm: ", trim(errmsg)
        error stop 1
      end if
    end if
    if (.not.self%config%is_configured) then
      log_fatal(*) "mLM configuration not set."
      error stop 1
    end if
    status = self%config%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "mLM config not valid: ", trim(errmsg)
      error stop 1
    end if

    self%output_active = .true.
    if (present(out_file)) then
      path = self%exchange%get_path(out_file, root=.true.)
      log_info(*) "Read mLM output config: ", path
      status = self%output_config%from_file(file=path, errmsg=errmsg)
      if (status /= NML_OK) then
        self%output_active = .false.
        log_warn(*) "mLM output disabled, config not found: ", trim(errmsg)
      end if
    end if
    if (self%output_config%is_configured) then
      status = self%output_config%is_valid(errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mLM output config invalid: ", trim(errmsg)
        error stop 1
      end if
    else
      self%output_active = .false.
      log_warn(*) "mLM output disabled, config not set."
    end if
    if (self%output_active) then
      status = self%config%is_set("output_path", idx=id, errmsg=errmsg)
      self%output_active = status == NML_OK
      if (status /= NML_OK) then
        log_warn(*) "mLM output disabled, path not set for domain ", n2s(id(1)), ": ", trim(errmsg)
      else
        self%output_path = self%exchange%get_path(self%config%output_path(id(1)))
      end if
    end if

    self%read_restart = self%config%read_restart(id(1))
    self%write_restart = self%config%write_restart(id(1))
    if (self%read_restart) then
      if (self%config%is_set("restart_input_path", idx=id) /= NML_OK) then
        log_fatal(*) "mLM: restart input path required when read_restart is true."
        error stop 1
      end if
      self%restart_input_path = self%exchange%get_path(self%config%restart_input_path(id(1)))
    end if
    if (self%write_restart) then
      if (self%config%is_set("restart_output_path", idx=id) /= NML_OK) then
        log_fatal(*) "mLM: restart output path required when write_restart is true."
        error stop 1
      end if
      self%restart_output_path = self%exchange%get_path(self%config%restart_output_path(id(1)))
    end if
  end subroutine mlm_configure

  !> \brief Resolve static lake metadata and publish hourly outflow.
  subroutine mlm_connect(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: n_lakes
    integer(i4) :: lake_case
    logical :: has_points, has_ids, has_levels, has_metadata

    lake_case = self%exchange%config%processes%lake
    has_points = associated(self%exchange%lake_points)
    has_ids = self%exchange%lake_ids%provided
    has_levels = self%exchange%lake_max_levels%provided
    has_metadata = has_points .and. has_ids .and. has_levels
    if ((has_points .neqv. has_ids) .or. (has_points .neqv. has_levels)) then
      log_fatal(*) "mLM: lake points, stable IDs, and maximum levels must be published together."
      error stop 1
    end if

    if (self%read_restart) call self%read_restart_metadata()
    if (has_metadata) then
      n_lakes = self%exchange%lake_points%n_points
      call self%exchange%lake_ids%require("mLM", .true., [n_lakes])
      call self%exchange%lake_max_levels%require("mLM", .true., [n_lakes])
      if (self%read_restart) call self%validate_restart_metadata()
      self%lake_ids = self%exchange%lake_ids%data
      call self%static_lake_points%init( &
        self%exchange%lake_points%x, self%exchange%lake_points%y, coordsys=self%exchange%lake_points%coordsys)
      self%static_lake_max_levels = self%exchange%lake_max_levels%data
    else
      if (.not.self%read_restart) then
        log_fatal(*) "mLM: lake process requires lake points."
        error stop 1
      end if
      n_lakes = self%restart_lake_points%n_points
      self%lake_ids = self%restart_lake_ids
      call self%static_lake_points%init( &
        self%restart_lake_points%x, self%restart_lake_points%y, coordsys=self%restart_lake_points%coordsys)
      self%static_lake_max_levels = self%restart_lake_max_levels
      self%exchange%lake_points => self%static_lake_points
      call self%exchange%lake_ids%publish_local("mLM", self%lake_ids, no_time)
      call self%exchange%lake_max_levels%publish_local("mLM", self%static_lake_max_levels, no_time)
      self%owns_restart_lake_points = .true.
      self%owns_restart_lake_metadata = .true.
    end if
    if (n_lakes < 1_i8) then
      log_fatal(*) "mLM: lake process requires at least one lake."
      error stop 1
    end if
    self%lake_ids = self%exchange%lake_ids%data
    if (lake_case == -2_i4) then
      if (associated(self%exchange%level0_lake) .and. self%exchange%lake_map%provided) then
        call self%exchange%lake_map%require("mLM", .true., [self%exchange%level0_lake%ncells])
        self%static_lake_grid = self%exchange%level0_lake
        self%lake_map = self%exchange%lake_map%data
      else if (self%read_restart) then
        call self%read_restart_topology()
        self%exchange%level0_lake => self%static_lake_grid
        call self%exchange%lake_map%publish_local("mLM", self%lake_map, no_time)
        self%owns_restart_lake_grid = .true.
      else
        log_fatal(*) "mLM process -2 requires a level-0 lake grid and lake map."
        error stop 1
      end if
      call lake_derive_area(self%static_lake_grid, self%lake_map, self%lake_ids, self%lake_area)
      call self%exchange%lake_area%publish_local("mLM", self%lake_area, no_time)
      call self%exchange%lake_pre%require("mLM", .true., [n_lakes])
      call self%exchange%lake_pet%require("mLM", .true., [n_lakes])
    end if
    allocate(self%outflow(n_lakes), source=0.0_dp)
    call self%exchange%lake_outflow%publish_local("mLM", self%outflow, 1_i4)
  end subroutine mlm_connect

  !> \brief Restore or initialize pass-through outflow after mRM has published inflow.
  subroutine mlm_initialize(self)
    class(mlm_t), target, intent(inout) :: self
    if (self%exchange%step_hours /= 1_i4) then
      log_fatal(*) "mLM lake processes require a one-hour model step."
      error stop 1
    end if
    call self%exchange%lake_inflow%require("mLM", .true., [size(self%lake_ids, kind=i8)])
    if (self%exchange%lake_inflow%stepping /= 1_i4) then
      log_fatal(*) "mLM: lake inflow must have one-hour support."
      error stop 1
    end if
    if (self%exchange%config%processes%lake == -2_i4) then
      call self%exchange%lake_area%require("mLM", .true., [size(self%lake_ids, kind=i8)])
    end if
    self%outflow = 0.0_dp
    if (self%read_restart) call self%read_restart_state()
    call self%validate_output_timing()
    call self%create_output()
  end subroutine mlm_initialize

  !> \brief Publish the preceding completed hourly inflow as current outflow.
  subroutine mlm_update(self)
    class(mlm_t), target, intent(inout) :: self
    self%outflow = self%exchange%lake_inflow%data
    if (self%exchange%config%processes%lake == -2_i4) then
      self%outflow = lake_balance_outflow(self%outflow, self%lake_area, self%exchange%lake_pre%data, &
        self%exchange%lake_pet%data, 3600.0_dp)
    end if
    call self%update_output()
  end subroutine mlm_update

  !> \brief Write optional restart and clear mLM-owned publication.
  subroutine mlm_finalize(self)
    class(mlm_t), target, intent(inout) :: self
    if (self%write_restart) call self%create_restart()
    if (self%output_active) then
      call self%ds_out%close()
      log_info(*) "Close mLM output file: ", self%output_path
    else
      log_info(*) "No mLM output file will be written"
    end if
    call self%exchange%lake_outflow%clear(owned=.true.)
    call self%exchange%lake_area%clear(owned=allocated(self%lake_area))
    if (self%owns_restart_lake_metadata) then
      call self%exchange%lake_ids%clear(owned=.true.)
      call self%exchange%lake_max_levels%clear(owned=.true.)
      self%owns_restart_lake_metadata = .false.
    end if
    if (self%owns_restart_lake_points) then
      if (associated(self%exchange%lake_points, self%static_lake_points)) nullify(self%exchange%lake_points)
      self%owns_restart_lake_points = .false.
    end if
    if (self%owns_restart_lake_grid) then
      call self%exchange%lake_map%clear(owned=.true.)
      if (associated(self%exchange%level0_lake, self%static_lake_grid)) nullify(self%exchange%level0_lake)
      self%owns_restart_lake_grid = .false.
    end if
    if (allocated(self%lake_ids)) deallocate(self%lake_ids)
    if (allocated(self%lake_map)) deallocate(self%lake_map)
    if (allocated(self%lake_area)) deallocate(self%lake_area)
    self%static_lake_points = points_t()
    if (allocated(self%static_lake_max_levels)) deallocate(self%static_lake_max_levels)
    self%restart_lake_points = points_t()
    if (allocated(self%restart_lake_ids)) deallocate(self%restart_lake_ids)
    if (allocated(self%restart_lake_max_levels)) deallocate(self%restart_lake_max_levels)
    if (allocated(self%outflow)) deallocate(self%outflow)
  end subroutine mlm_finalize

  !> \brief Create the lake-point output dataset in stable exchange order.
  subroutine mlm_create_output(self)
    class(mlm_t), target, intent(inout) :: self
    type(var), allocatable :: vars(:)
    character(:), allocatable :: delta, dtype
    integer(i4) :: timestamp

    if (.not.self%output_active) return

    timestamp = self%output_config%output_time_reference
    delta = time_units_delta(self%output_config%output_frequency, timestamp)
    dtype = "f64"
    if (.not.self%output_config%output_double_precision) dtype = "f32"

    vars = [var(name="lake_id", long_name="stable lake ID", dtype="i64", kind="i8", static=.true.)]
    if (self%output_config%out_lake_outflow) then
      vars = [vars, self%exchange%lake_outflow%as_output_var(dtype=dtype, avg=.true.)]
    end if

    log_info(*) "Create mLM lake-point output file: ", self%output_path
    call self%ds_out%init( &
      path          = self%output_path, &
      points        = self%exchange%lake_points, &
      vars          = vars, &
      start_time    = self%exchange%start_time, &
      delta         = delta, &
      timestamp     = timestamp, &
      deflate_level = self%output_config%output_deflate_level, &
      point_dim_name = "lake", &
      time_series   = .true.)
    call self%ds_out%update("lake_id", self%lake_ids)
    call self%ds_out%write_static()
  end subroutine mlm_create_output

  !> \brief Validate mLM output cadence and restart/output boundary alignment.
  subroutine mlm_validate_output_timing(self)
    class(mlm_t), target, intent(in) :: self
    integer(i4) :: frequency

    if (.not.self%output_active) return
    frequency = self%output_config%output_frequency
    select case (frequency)
      case (daily, monthly, yearly, no_time)
        continue
      case default
        if (frequency < 1_i4 .or. mod(frequency, self%exchange%step_hours) /= 0_i4) then
          log_fatal(*) "mLM output_frequency=", frequency, &
            "h must be a positive whole multiple of model_step=", self%exchange%step_hours, "h."
          error stop 1
        end if
    end select

    if (self%write_restart .and. .not.self%at_output_boundary(self%exchange%end_time)) then
      log_fatal(*) "mLM restart output time must coincide with an mLM output boundary."
      error stop 1
    end if
  end subroutine mlm_validate_output_timing

  !> \brief Report whether a timestamp is a configured mLM output boundary.
  logical function mlm_at_output_boundary(self, time) result(boundary)
    class(mlm_t), target, intent(in) :: self
    type(datetime), intent(in) :: time
    integer(i4) :: elapsed_hours

    select case (self%output_config%output_frequency)
      case (daily)
        boundary = time%is_new_day()
      case (monthly)
        boundary = time%is_new_month()
      case (yearly)
        boundary = time%is_new_year()
      case (no_time)
        boundary = time == self%exchange%end_time
      case default
        elapsed_hours = nint((time - self%exchange%start_time) / one_hour(), i4)
        boundary = elapsed_hours >= 0_i4 .and. &
          mod(elapsed_hours, self%output_config%output_frequency) == 0_i4
    end select
  end function mlm_at_output_boundary

  !> \brief Buffer the current lake release and write completed output intervals.
  subroutine mlm_update_output(self)
    class(mlm_t), target, intent(inout) :: self

    if (.not.self%output_active) return
    if (self%output_config%out_lake_outflow) then
      call self%ds_out%update("lake_outflow", self%outflow)
    end if
    if (self%at_output_boundary(self%exchange%time)) call self%ds_out%write(self%exchange%time)
  end subroutine mlm_update_output

  !> \brief Persist static lake metadata and current outflow with stable IDs.
  subroutine mlm_create_restart(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcDimension) :: dims(0), lake_dim
    type(NcVariable) :: var

    if (self%static_lake_points%n_points /= size(self%lake_ids, kind=i8) .or. &
        .not.allocated(self%static_lake_max_levels)) then
      log_fatal(*) "mLM: cannot write restart without static lake metadata."
      error stop 1
    end if
    nc = NcDataset(self%restart_output_path, "w")
    lake_dim = nc%setDimension("lake", int(size(self%lake_ids), i4))
    var = nc%setVariable("lake_id", "i64", [lake_dim])
    call var%setData(self%lake_ids)
    var = nc%setVariable("lake_outlet_x", "f64", [lake_dim])
    call var%setData(self%static_lake_points%x)
    var = nc%setVariable("lake_outlet_y", "f64", [lake_dim])
    call var%setData(self%static_lake_points%y)
    var = nc%setVariable("lake_max_level", "f64", [lake_dim])
    call var%setAttribute("units", "m")
    call var%setData(self%static_lake_max_levels)
    var = nc%setVariable("lake_coordsys", "i32", dims(:0))
    call var%setData(self%static_lake_points%coordsys)
    var = nc%setVariable("lake_outflow", "f64", [lake_dim])
    call var%setAttribute("units", "m3 s-1")
    call var%setData(self%outflow)
    if (self%exchange%config%processes%lake == -2_i4) then
      call lake_write_topology(nc, self%static_lake_grid, self%lake_map)
    end if
    var = nc%setVariable("mlm_meta", "i32", dims(:0))
    call var%setData(0_i4)
    call var%setAttribute("time_stamp", self%exchange%time%str())
    call nc%close()
  end subroutine mlm_create_restart

  !> \brief Read and validate static lake metadata before resolving the exchange provider.
  subroutine mlm_read_restart_metadata(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcVariable) :: var, meta_var
    real(dp), allocatable :: outlet_x(:), outlet_y(:)
    integer(i4), allocatable :: coordsys_data
    integer(i8) :: i
    character(64) :: restart_time
    character(:), allocatable :: expected_time

    nc = NcDataset(self%restart_input_path, "r")
    if (.not.nc%hasVariable("lake_id") .or. .not.nc%hasVariable("lake_outflow")) then
      log_fatal(*) "mLM restart is missing lake IDs or lake outflow: ", self%restart_input_path
      error stop 1
    end if
    var = nc%getVariable("lake_id")
    call var%getData(self%restart_lake_ids)
    if (size(self%restart_lake_ids, kind=i8) < 1_i8) then
      log_fatal(*) "mLM restart contains no lakes: ", self%restart_input_path
      error stop 1
    end if
    if (any(self%restart_lake_ids <= 0_i8)) then
      log_fatal(*) "mLM restart contains a non-positive stable lake ID."
      error stop 1
    end if
    do i = 1_i8, size(self%restart_lake_ids, kind=i8) - 1_i8
      if (any(self%restart_lake_ids(i + 1_i8:) == self%restart_lake_ids(i))) then
        log_fatal(*) "mLM restart contains duplicate stable lake IDs."
        error stop 1
      end if
    end do

    if (.not.nc%hasVariable("lake_outlet_x") .or. .not.nc%hasVariable("lake_outlet_y") .or. &
        .not.nc%hasVariable("lake_max_level") .or. .not.nc%hasVariable("lake_coordsys")) then
      log_fatal(*) "mLM restart is missing required static lake metadata: ", self%restart_input_path
      error stop 1
    end if
    var = nc%getVariable("lake_outlet_x")
    call var%getData(outlet_x)
    var = nc%getVariable("lake_outlet_y")
    call var%getData(outlet_y)
    var = nc%getVariable("lake_max_level")
    call var%getData(self%restart_lake_max_levels)
    var = nc%getVariable("lake_coordsys")
    call var%getData(coordsys_data)
    if (size(outlet_x, kind=i8) /= size(self%restart_lake_ids, kind=i8) .or. &
        size(outlet_y, kind=i8) /= size(self%restart_lake_ids, kind=i8) .or. &
        size(self%restart_lake_max_levels, kind=i8) /= size(self%restart_lake_ids, kind=i8)) then
      log_fatal(*) "mLM restart static lake metadata does not align with stable lake IDs."
      error stop 1
    end if
    if (coordsys_data /= cartesian .and. coordsys_data /= spherical) then
      log_fatal(*) "mLM restart contains an unsupported lake coordinate system."
      error stop 1
    end if
    if (.not.all(ieee_is_finite(outlet_x)) .or. .not.all(ieee_is_finite(outlet_y)) .or. &
        .not.all(ieee_is_finite(self%restart_lake_max_levels))) then
      log_fatal(*) "mLM restart static lake metadata must be finite."
      error stop 1
    end if
    call self%restart_lake_points%init(outlet_x, outlet_y, coordsys=coordsys_data)

    if (.not.nc%hasVariable("mlm_meta")) then
      log_fatal(*) "mLM restart metadata variable mlm_meta is missing: ", self%restart_input_path
      error stop 1
    end if
    meta_var = nc%getVariable("mlm_meta")
    if (.not.meta_var%hasAttribute("time_stamp")) then
      log_fatal(*) "mLM restart metadata has no time_stamp: ", self%restart_input_path
      error stop 1
    end if
    call meta_var%getAttribute("time_stamp", restart_time)
    expected_time = self%exchange%start_time%str()
    if (trim(restart_time) /= trim(expected_time)) then
      log_fatal(*) "mLM restart timestamp ", trim(restart_time), " does not match domain restart time ", trim(expected_time), "."
      error stop 1
    end if
    call nc%close()
  end subroutine mlm_read_restart_metadata

  !> \brief Restore the static level-0 lake footprint required by lake case -2.
  subroutine mlm_read_restart_topology(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    nc = NcDataset(self%restart_input_path, "r")
    call lake_read_topology(nc, self%static_lake_grid, self%lake_map)
    call nc%close()
  end subroutine mlm_read_restart_topology

  !> \brief Verify restart metadata against Input metadata without imposing restart-file ordering.
  subroutine mlm_validate_restart_metadata(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: i, j, n_lakes

    n_lakes = self%exchange%lake_points%n_points
    if (size(self%restart_lake_ids, kind=i8) /= n_lakes) then
      log_fatal(*) "mLM restart lake IDs do not match Input lake metadata."
      error stop 1
    end if
    if (self%restart_lake_points%coordsys /= self%exchange%lake_points%coordsys) then
      log_fatal(*) "mLM restart and Input lake coordinate systems differ."
      error stop 1
    end if
    do i = 1_i8, n_lakes - 1_i8
      if (any(self%exchange%lake_ids%data(i + 1_i8:) == self%exchange%lake_ids%data(i))) then
        log_fatal(*) "Input lake metadata contains duplicate stable lake IDs."
        error stop 1
      end if
    end do
    do i = 1_i8, n_lakes
      do j = 1_i8, size(self%restart_lake_ids, kind=i8)
        if (self%restart_lake_ids(j) == self%exchange%lake_ids%data(i)) exit
      end do
      if (j > size(self%restart_lake_ids, kind=i8)) then
        log_fatal(*) "mLM restart and Input lake metadata have different stable lake IDs."
        error stop 1
      end if
      if (.not.is_close(self%restart_lake_points%x(j), self%exchange%lake_points%x(i)) .or. &
          .not.is_close(self%restart_lake_points%y(j), self%exchange%lake_points%y(i)) .or. &
          .not.is_close(self%restart_lake_max_levels(j), self%exchange%lake_max_levels%data(i))) then
        log_fatal(*) "mLM restart and Input lake metadata differ for stable lake ID."
        error stop 1
      end if
    end do
  end subroutine mlm_validate_restart_metadata

  !> \brief Restore outflow by stable ID rather than restart-file order.
  subroutine mlm_read_restart_state(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcVariable) :: var
    integer(i8), allocatable :: ids(:)
    real(dp), allocatable :: values(:)
    integer(i8) :: i, j
    nc = NcDataset(self%restart_input_path, "r")
    if (.not.nc%hasVariable("lake_id") .or. .not.nc%hasVariable("lake_outflow")) then
      log_fatal(*) "mLM restart is missing lake IDs or lake outflow: ", self%restart_input_path
      error stop 1
    end if
    var = nc%getVariable("lake_id")
    call var%getData(ids)
    var = nc%getVariable("lake_outflow")
    call var%getData(values)
    call nc%close()
    if (size(ids, kind=i8) /= size(self%lake_ids, kind=i8) .or. size(values, kind=i8) /= size(ids, kind=i8)) then
      log_fatal(*) "mLM restart lake count does not match configured lakes."
      error stop 1
    end if
    do i = 1_i8, size(ids, kind=i8) - 1_i8
      if (any(ids(i + 1_i8:) == ids(i))) then
        log_fatal(*) "mLM restart contains duplicate stable lake IDs."
        error stop 1
      end if
    end do
    do i = 1_i8, size(self%lake_ids, kind=i8)
      do j = 1_i8, size(ids, kind=i8)
        if (ids(j) == self%lake_ids(i)) exit
      end do
      if (j > size(ids, kind=i8)) then
        log_fatal(*) "mLM restart is missing stable lake ID."
        error stop 1
      end if
      self%outflow(i) = values(j)
    end do
  end subroutine mlm_read_restart_state

end module mo_mlm_container
