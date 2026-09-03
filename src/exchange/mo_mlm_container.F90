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
  use mo_exchange_type, only: exchange_t
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable
  use mo_grid_io, only: no_time
  use mo_points, only: points_t, cartesian, spherical
  use mo_utils, only: is_close
  use mo_string_utils, only: n2s => num2str
  use nml_config_lake, only: nml_config_lake_t
  use nml_helper, only: NML_OK
  implicit none

  character(*), parameter :: s = "mlm" !< logging scope

  !> \class mlm_t
  !> \brief Domain-owned lake pass-through component.
  !> \authors Pallav Shrestha
  !> \authors Sebastian Mueller
  type, public :: mlm_t
    type(nml_config_lake_t) :: config !< lake configuration
    type(exchange_t), pointer :: exchange => null() !< owning domain exchange
    integer(i8), allocatable :: lake_ids(:) !< stable IDs in point-set order
    type(points_t) :: static_lake_points !< static point metadata retained for restart output
    real(dp), allocatable :: static_lake_max_levels(:) !< maximum levels retained for restart output
    type(points_t) :: restart_lake_points !< point set restored from restart metadata
    integer(i8), allocatable :: restart_lake_ids(:) !< stable IDs in restart-file order
    real(dp), allocatable :: restart_lake_max_levels(:) !< maximum levels restored from restart metadata
    real(dp), allocatable :: outflow(:) !< current lake outflow [m3 s-1]
    logical :: active = .false. !< whether process -1 is selected
    logical :: read_restart = .false. !< read mLM restart
    logical :: write_restart = .false. !< write mLM restart
    logical :: owns_restart_lake_points = .false. !< whether mLM published the restart point set
    logical :: owns_restart_lake_metadata = .false. !< whether mLM published restart IDs and maximum levels
    character(:), allocatable :: restart_input_path !< resolved input restart path
    character(:), allocatable :: restart_output_path !< resolved output restart path
  contains
    procedure :: set_dims => mlm_set_dims
    procedure :: configure => mlm_configure
    procedure :: connect => mlm_connect
    procedure :: initialize => mlm_initialize
    procedure :: update => mlm_update
    procedure :: finalize => mlm_finalize
    procedure, private :: create_restart => mlm_create_restart
    procedure, private :: read_restart_metadata => mlm_read_restart_metadata
    procedure, private :: read_restart_state => mlm_read_restart_state
    procedure, private :: validate_restart_metadata => mlm_validate_restart_metadata
  end type mlm_t

contains

  !> \brief Set generated lake-configuration dimensions.
  subroutine mlm_set_dims(self)
    class(mlm_t), target, intent(inout) :: self
    character(1024) :: errmsg
    integer :: status
    status = self%config%set_dims(n_domains=self%exchange%nml_n_domains, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "mLM: error setting config_lake dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine mlm_set_dims

  !> \brief Configure process -1 and optional restart paths.
  subroutine mlm_configure(self, file)
    class(mlm_t), target, intent(inout) :: self
    character(*), optional, intent(in) :: file
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer(i4) :: id(1), lake_case
    integer :: status

    lake_case = self%exchange%config%processes%lake
    self%active = lake_case /= 0_i4
    if (lake_case /= 0_i4 .and. lake_case /= -1_i4) then
      log_fatal(*) "mLM: unsupported lake process case: ", n2s(lake_case)
      error stop 1
    end if
    if (.not.self%active) return
    if (self%exchange%config%processes%routing == 0_i4) then
      log_fatal(*) "mLM process -1 requires active mRM routing."
      error stop 1
    end if
    id(1) = self%exchange%nml_domain_id
    if (present(file)) then
      path = self%exchange%get_path(file)
      status = self%config%from_file(path, errmsg=errmsg)
      if (status /= NML_OK) then
        log_fatal(*) "mLM: error reading config_lake: ", trim(errmsg)
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
    logical :: has_points, has_ids, has_levels, has_metadata

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
    allocate(self%outflow(n_lakes), source=0.0_dp)
    call self%exchange%lake_outflow%publish_local("mLM", self%outflow, 1_i4)
  end subroutine mlm_connect

  !> \brief Restore or initialize pass-through outflow after mRM has published inflow.
  subroutine mlm_initialize(self)
    class(mlm_t), target, intent(inout) :: self
    if (self%exchange%step_hours /= 1_i4) then
      log_fatal(*) "mLM process -1 requires a one-hour model step."
      error stop 1
    end if
    call self%exchange%lake_inflow%require("mLM", .true., [size(self%lake_ids, kind=i8)])
    if (self%exchange%lake_inflow%stepping /= 1_i4) then
      log_fatal(*) "mLM: lake inflow must have one-hour support."
      error stop 1
    end if
    if (self%read_restart) call self%read_restart_state()
  end subroutine mlm_initialize

  !> \brief Publish the preceding completed hourly inflow as current outflow.
  subroutine mlm_update(self)
    class(mlm_t), target, intent(inout) :: self
    self%outflow = self%exchange%lake_inflow%data
  end subroutine mlm_update

  !> \brief Write optional restart and clear mLM-owned publication.
  subroutine mlm_finalize(self)
    class(mlm_t), target, intent(inout) :: self
    if (self%write_restart) call self%create_restart()
    call self%exchange%lake_outflow%clear(owned=.true.)
    if (self%owns_restart_lake_metadata) then
      call self%exchange%lake_ids%clear(owned=.true.)
      call self%exchange%lake_max_levels%clear(owned=.true.)
      self%owns_restart_lake_metadata = .false.
    end if
    if (self%owns_restart_lake_points) then
      if (associated(self%exchange%lake_points, self%static_lake_points)) nullify(self%exchange%lake_points)
      self%owns_restart_lake_points = .false.
    end if
    if (allocated(self%lake_ids)) deallocate(self%lake_ids)
    self%static_lake_points = points_t()
    if (allocated(self%static_lake_max_levels)) deallocate(self%static_lake_max_levels)
    self%restart_lake_points = points_t()
    if (allocated(self%restart_lake_ids)) deallocate(self%restart_lake_ids)
    if (allocated(self%restart_lake_max_levels)) deallocate(self%restart_lake_max_levels)
    if (allocated(self%outflow)) deallocate(self%outflow)
  end subroutine mlm_finalize

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
    var = nc%getVariable("lake_outlet_x"); call var%getData(outlet_x)
    var = nc%getVariable("lake_outlet_y"); call var%getData(outlet_y)
    var = nc%getVariable("lake_max_level"); call var%getData(self%restart_lake_max_levels)
    var = nc%getVariable("lake_coordsys"); call var%getData(coordsys_data)
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
    var = nc%getVariable("lake_id"); call var%getData(ids)
    var = nc%getVariable("lake_outflow"); call var%getData(values)
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
