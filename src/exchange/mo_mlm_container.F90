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
  use mo_logging
  use mo_kind, only: i4, i8, dp
  use mo_exchange_type, only: exchange_t
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable
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
    real(dp), allocatable :: outflow(:) !< current lake outflow [m3 s-1]
    logical :: active = .false. !< whether process -1 is selected
    logical :: read_restart = .false. !< read mLM restart
    logical :: write_restart = .false. !< write mLM restart
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
    procedure, private :: read_restart_state => mlm_read_restart_state
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

  !> \brief Validate static lake metadata and publish hourly outflow.
  subroutine mlm_connect(self)
    class(mlm_t), target, intent(inout) :: self
    integer(i8) :: n_lakes
    if (.not.associated(self%exchange%lake_points)) then
      log_fatal(*) "mLM: lake process requires lake points."
      error stop 1
    end if
    n_lakes = self%exchange%lake_points%n_points
    call self%exchange%lake_ids%require("mLM", .true., [n_lakes])
    call self%exchange%lake_max_levels%require("mLM", .true., [n_lakes])
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
    if (allocated(self%lake_ids)) deallocate(self%lake_ids)
    if (allocated(self%outflow)) deallocate(self%outflow)
  end subroutine mlm_finalize

  !> \brief Persist current outflow with stable IDs.
  subroutine mlm_create_restart(self)
    class(mlm_t), target, intent(inout) :: self
    type(NcDataset) :: nc
    type(NcDimension) :: lake_dim
    type(NcVariable) :: var
    nc = NcDataset(self%restart_output_path, "w")
    lake_dim = nc%setDimension("lake", int(size(self%lake_ids), i4))
    var = nc%setVariable("lake_id", "i64", [lake_dim])
    call var%setData(self%lake_ids)
    var = nc%setVariable("lake_outflow", "f64", [lake_dim])
    call var%setAttribute("units", "m3 s-1")
    call var%setData(self%outflow)
    call nc%close()
  end subroutine mlm_create_restart

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
