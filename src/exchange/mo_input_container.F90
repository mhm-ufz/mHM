!> \file    mo_input_container.f90
!> \brief   \copybrief mo_input_container
!> \details \copydetails mo_input_container

!> \brief   Module for an input container.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Aug 2025
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_exchange
#include "logging.h"
module mo_input_container
  use mo_logging
  use mo_kind, only: i2, i4, i8, dp
  use mo_list, only: list
  use mo_os, only: path_ext
  use mo_exchange_type, only: exchange_t, var_dp
  use mo_datetime, only: datetime, timedelta, HOUR_SECONDS, DAY_HOURS, one_hour, one_day
  use mo_grid, only: grid_t, data_t, cartesian, spherical
  use mo_river, only: river_t
  use mo_points, only: points_t
  use mo_points_io, only: points_input_dataset
  use mo_grid_io, only: var, input_dataset, end_timestamp, start_timestamp, no_time, daily, monthly, yearly, varying
  use mo_string_utils, only: n2s => num2str
  use nml_config_input, only: nml_config_input_t
  use nml_config_coupling, only: nml_config_coupling_t
  use nml_helper, only: NML_OK

  character(len=*), parameter :: s = "input" !< module scope for logging

  !> \class   input_var_abc
  !> \brief   Abstract Base Class for a single input variable.
  type, abstract :: input_var_abc
    logical :: provided = .false. !< whether this input variable is provided
    logical :: coupled = .false. !< whether this variable is from a coupler
    logical :: static = .false. !< whether this variable is static in time
    logical :: allow_static = .false. !< whether a static file is accepted for a temporal variable
    logical :: morph_latlon = .false. !< whether morphology grid should be initialized as spherical coordinates
    character(len=:), allocatable :: path !< path to the input dataset
    character(len=:), allocatable :: name !< variable name in the dataset
    type(input_dataset) :: ds !< input netcdf dataset
    type(grid_t), pointer :: grid => null() !< grid reference used for data reads
    integer(i4) :: var_id !< variable ID in the dataset
    integer(i4) :: stepping = 0_i4 !< time-stepping in hours (0 - static, -1 - daily, -2 -monthly, -3 - yearly, >0 - hourly)
    integer(i4) :: offset !< offset for reading in chunked mode
    type(datetime) :: chunk_time_start !< start time of current chunk (time span start when coupled)
    type(datetime) :: chunk_time_end !< end time of current chunk (time span end when coupled)
    logical, allocatable :: mask(:,:) !< optional mask for this variable size(nx, ny) (used for coupling)
  contains
    procedure :: init => input_var_init
    procedure :: is_ascii => input_var_is_ascii
    procedure :: open_dataset => input_var_open_dataset
    procedure :: check_couple_status => input_var_check_couple_status
    procedure :: set_mask => input_var_set_mask
  end type input_var_abc

  !> \class   input_var_dp
  !> \brief   Class for a single double precision input variable.
  type, public, extends(input_var_abc) :: input_var_dp
    real(dp), allocatable :: cache(:,:) !< data array for this variable size(ncells, ntimes)
  contains
    procedure :: update => input_var_dp_update
    procedure :: read_static => input_var_dp_read_static
    procedure :: reset_time => input_var_dp_reset_time
    ! procedure :: set => input_var_dp_set
  end type input_var_dp

  !> \class   input_var_i4
  !> \brief   Class for a single integer input variable.
  type, public, extends(input_var_abc) :: input_var_i4
    integer(i4), allocatable :: cache(:,:) !< data array for this variable size(ncells, ntimes)
  contains
    procedure :: update => input_var_i4_update
    procedure :: read_static => input_var_i4_read_static
    procedure :: reset_time => input_var_i4_reset_time
    ! procedure :: set => input_var_i4_set
  end type input_var_i4

  !> \class   input_var_i2
  !> \brief   Class for a single 16-bit integer input variable.
  type, public, extends(input_var_abc) :: input_var_i2
    integer(i2), allocatable :: cache(:,:) !< data array for this variable size(ncells, ntimes)
  contains
    procedure :: update => input_var_i2_update
    procedure :: read_static => input_var_i2_read_static
    procedure :: reset_time => input_var_i2_reset_time
  end type input_var_i2

  !> \class   input_var2d_dp
  !> \brief   Class for a single 2D double precision input variable.
  type, public, extends(input_var_abc) :: input_var2d_dp
    real(dp), allocatable :: cache(:,:,:) !< data array for this variable size(ncells, nlayers, ntimes)
  contains
    procedure :: update => input_var2d_dp_update
    procedure :: read_static => input_var2d_dp_read_static
    procedure :: reset_time => input_var2d_dp_reset_time
    ! procedure :: set => input_var2d_dp_set
  end type input_var2d_dp

  !> \class   input_var2d_i4
  !> \brief   Class for a single 2D integer input variable.
  type, public, extends(input_var_abc) :: input_var2d_i4
    integer(i4), allocatable :: cache(:,:,:) !< data array for this variable size(ncells, nlayers, ntimes)
  contains
    procedure :: update => input_var2d_i4_update
    procedure :: read_static => input_var2d_i4_read_static
    procedure :: reset_time => input_var2d_i4_reset_time
    ! procedure :: set => input_var2d_i4_set
  end type input_var2d_i4

  !> \class   input_config_t
  !> \brief   Class for a single Input container.
  type, public :: input_config_t
    logical :: active = .true. !< flag to activate the Input container
    integer(i4) :: domain !< domain number to read correct configuration
   type(nml_config_input_t) :: input !< configuration for inputs from files
   type(nml_config_coupling_t) :: coupling !< configuration for coupling
  contains
    procedure :: read => input_config_read
  end type input_config_t

  !> \class   input_t
  !> \brief   Class for a single Input container.
  type, public :: input_t
    type(input_config_t) :: config !< configuration of the Input container
    type(exchange_t), pointer :: exchange => null() !< exchange container of the domain
    type(grid_t) :: tgt_level0 !< grid level 0 of the domain if given from input
    type(grid_t) :: tgt_level0_land !< level-0 land grid derived by excluding lake cells
    type(grid_t) :: tgt_level0_lake !< level-0 lake grid derived from lake footprints
    type(grid_t) :: tgt_level1 !< grid level 1 of the domain if given from input
    type(grid_t) :: tgt_level2 !< grid level 2 of the domain if given from input
    type(grid_t) :: tgt_level3 !< grid level 3 of the domain if given from input
    type(river_t) :: river_l0 !< full level-0 river network derived from flow direction
    logical :: owns_river_l0 = .false. !< whether this input instance published the level-0 river
    type(points_t) :: lake_outlets !< configured lake outlet coordinates
    integer(i8), allocatable :: lake_ids(:) !< stable lake IDs aligned with lake_outlets
    real(dp), allocatable :: lake_max_levels(:) !< maximum lake levels aligned with lake_outlets
    integer(i4) :: chunking !< chunking configuration (0 single read, -1 daily, -2 monthly, -3 yearly, >0 every n hours)
    integer(i4) :: time_stamp_location !< location of time-stamp variable in input datasets (0 start, 1 center, 2 end)
    logical :: morph_latlon = .false. !< whether morphology inputs are defined in spherical (lat/lon) coordinates
    ! level 0 inputs
    type(input_var_i4) :: morph_mask !< input variable for morph mask
    type(input_var_dp) :: dem !< input variable for digital elevation model
    type(input_var_dp) :: slope !< input variable for slope
    type(input_var_dp) :: aspect !< input variable for aspect
    type(input_var_i2) :: fdir !< input variable for flow direction
    type(input_var_i4) :: facc !< input variable for flow accumulation
    type(input_var_i4) :: geo_class !< input variable for geology class
    type(input_var_i4) :: soil_class !< input variable for soil class
    type(input_var2d_i4) :: soil_horizon_class !< input variable for soil horizon class
    type(input_var_i4) :: lai_class !< input variable for LAI class
    type(input_var_i4) :: meteo_mask !< input variable for the meteorological mask
    type(input_var_dp) :: pre !< raw precipitation input variable on level2
    type(input_var_dp) :: pet !< raw PET input variable on level2
    type(input_var_dp) :: temp !< raw temperature input variable on level2
    type(input_var_dp) :: tann !< raw annual mean temperature input variable on level2
    type(input_var_dp) :: tmin !< raw minimum temperature input variable on level2
    type(input_var_dp) :: tmax !< raw maximum temperature input variable on level2
    type(input_var_dp) :: ssrd !< raw short-wave radiation input variable on level2
    type(input_var_dp) :: strd !< raw long-wave radiation input variable on level2
    type(input_var_dp) :: netrad !< raw net radiation input variable on level2
    type(input_var_dp) :: eabs !< raw vapor pressure input variable on level2
    type(input_var_dp) :: wind !< raw wind input variable on level2
    ! level 1 inputs
    type(input_var_i4) :: hydro_mask !< input variable for hydro mask
    type(input_var_dp) :: runoff !< input variable for runoff
    integer(i4), allocatable :: soil_class_one_layer(:,:) !< adapter for single-layer soil class input
  contains
    procedure :: set_dims => input_set_dims
    procedure :: configure => input_configure
    procedure :: connect => input_connect
    procedure :: initialize => input_initialize
    procedure :: update => input_update
    procedure :: finalize => input_finalize
    procedure, private :: read_fdir_file => input_read_fdir_file
    procedure, private :: build_river_l0 => input_build_river_l0
    procedure, private :: read_lake_specification => input_read_lake_specification
    procedure, private :: build_lake_grids => input_build_lake_grids
  end type input_t

contains

  !> \brief Check whether the grid needs to be initialized in the exchange variable.
  logical function need_grid(tgt_grid, exchange_grid)
    type(grid_t), intent(in), target :: tgt_grid !< target grid in input
    type(grid_t), intent(inout), pointer :: exchange_grid !< grid pointer in exchange
    need_grid = .not.associated(exchange_grid)
    if (need_grid) exchange_grid => tgt_grid
  end function need_grid

  !> \brief Connect the land grid, aliasing the full grid when no derived lake mask exists.
  logical function need_level0_land_grid(tgt_grid, full_grid, land_grid)
    type(grid_t), intent(in), target :: tgt_grid
    type(grid_t), intent(inout), pointer :: full_grid
    type(grid_t), intent(inout), pointer :: land_grid
    need_level0_land_grid = .false.
    if (associated(land_grid)) return
    need_level0_land_grid = need_grid(tgt_grid, full_grid)
    land_grid => full_grid
  end function need_level0_land_grid

  !> \brief Repack a full level-0 real cache onto the land grid without rereading its source.
  subroutine repack_l0_dp_cache(full_grid, land_grid, cache, name)
    type(grid_t), intent(in) :: full_grid
    type(grid_t), intent(in) :: land_grid
    real(dp), allocatable, intent(inout) :: cache(:,:)
    character(*), intent(in) :: name
    real(dp), allocatable :: data_2d(:,:), land_cache(:,:)
    integer(i4) :: layer

    if (.not.allocated(cache)) then
      log_fatal(*) "Input: cannot repack uninitialized level-0 data: ", name
      error stop 1
    end if
    if (size(cache, 1, kind=i8) /= full_grid%ncells) then
      log_fatal(*) "Input: full-grid data size does not match level-0 cells before land repacking: ", name
      error stop 1
    end if
    allocate(data_2d(full_grid%nx, full_grid%ny))
    allocate(land_cache(land_grid%ncells, size(cache, 2)))
    do layer = 1_i4, size(cache, 2)
      call full_grid%unpack_into(cache(:, layer), data_2d)
      call land_grid%pack_into(data_2d, land_cache(:, layer))
    end do
    call move_alloc(land_cache, cache)
  end subroutine repack_l0_dp_cache

  !> \brief Copy stepping/static/provided metadata from an input variable to the exchange variable.
  subroutine sync_input_var_meta(input_var, exchange_var)
    class(input_var_abc), intent(in) :: input_var
    type(var_dp), intent(inout) :: exchange_var
    exchange_var%provided = input_var%provided
    exchange_var%static = input_var%static
    call exchange_var%set_stepping("Input", input_var%stepping)
  end subroutine sync_input_var_meta

  !> \brief Update the time frame for chunked reading.
  subroutine update_time_frame(chunk_time_start, chunk_time_end, chunking, end_time)
    type(datetime), intent(inout) :: chunk_time_start !< start time of current chunk
    type(datetime), intent(inout) :: chunk_time_end !< end time of current chunk
    integer(i4), intent(in) :: chunking !< chunking configuration (0 single read, -1 daily, -2 monthly, -3 yearly, >0 every n hours)
    type(datetime), intent(in) :: end_time !< end time of the simulation (for chunking limits)
    chunk_time_start = chunk_time_end
    select case(chunking)
      case(daily)
        chunk_time_end = chunk_time_start%next_new_day()
      case(monthly)
        chunk_time_end = chunk_time_start%next_new_month()
      case(yearly)
        chunk_time_end = chunk_time_start%next_new_year()
      case(0_i4) ! flag for reading all at once
        chunk_time_end = end_time ! read once until end time
      case default
        if (chunking < 1_i4) then
          log_fatal(*) "Input: Chunk not valid: ", n2s(chunking)
          error stop 1
        end if
        chunk_time_end = chunk_time_start + timedelta(hours=chunking)
      end select
    if (chunk_time_end > end_time) chunk_time_end = end_time
  end subroutine update_time_frame

  !> \brief Initialize the input variable.
  subroutine input_var_init(self, path, name, static, coupled, morph_latlon, allow_static)
    class(input_var_abc), intent(inout), target :: self
    character(*), intent(in), optional :: path !< path to the input dataset
    character(*), intent(in), optional :: name !< variable name in the dataset
    logical, intent(in), optional :: static !< whether this variable is static in time
    logical, intent(in), optional :: coupled !< whether this variable is from a coupler
    logical, intent(in), optional :: morph_latlon !< whether to initialize ASCII grid as spherical coordinates
    logical, intent(in), optional :: allow_static !< whether to allow a static input dataset for a temporal variable
    self%provided = .true. ! mark as provided when initialized
    if (present(path)) self%path = trim(path)
    if (present(name)) self%name = trim(name)
    if (present(static)) self%static = static
    if (present(coupled)) self%coupled = coupled
    if (present(morph_latlon)) self%morph_latlon = morph_latlon
    if (present(allow_static)) self%allow_static = allow_static
  end subroutine input_var_init

  logical function input_var_is_ascii(self)
    class(input_var_abc), intent(in), target :: self
    character(:), allocatable :: ext
    input_var_is_ascii = .false.
    if (.not.allocated(self%path)) return
    ext = path_ext(self%path)
    input_var_is_ascii = (ext == ".asc") .or. (ext == ".ASC")
  end function input_var_is_ascii

  subroutine input_var_open_dataset(self, kind, timestamp, grid, init_grid, layered)
    class(input_var_abc), intent(inout), target :: self
    character(*), intent(in) :: kind !< kind of variable to create dataset for
    integer(i4), intent(in) :: timestamp !< time stamp for the dataset
    type(grid_t), intent(inout), pointer :: grid !< grid for the variable
    logical, intent(in) :: init_grid !< whether to initialize the grid
    logical, intent(in), optional :: layered !< whether the variable has a layer dimension
    character(:), allocatable :: grid_init_var
    logical :: layered_
    if (.not.self%provided) return
    layered_ = .false.
    if (present(layered)) layered_ = layered
    self%grid => grid

    scope_info(s,*) "Prepare input variable '", trim(self%name), "' from file: ", trim(self%path)
    if (init_grid) then
      scope_info(s,*) "Initialize input grid from '", trim(self%name), "'"
    end if

    if (self%is_ascii()) then
      if (.not.self%static) then
        log_fatal(*) "Input: ASCII input is only supported for static variables. Variable: ", trim(self%name)
        error stop 1
      end if
      if (layered_) then
        log_fatal(*) "Input: layered ASCII input is not supported. Variable: ", trim(self%name)
        error stop 1
      end if
      if (init_grid) then
        call self%grid%from_ascii_file(self%path, coordsys=merge(spherical, cartesian, self%morph_latlon))
      else if (.not.associated(self%grid)) then
        log_fatal(*) "Input: grid is not connected for ASCII variable: ", trim(self%name)
        error stop 1
      end if
      return
    end if

    if (init_grid) then
      grid_init_var = self%name ! allocate variable name for grid initialization
      call self%ds%init(path=self%path, vars=[var(name=trim(self%name), kind=kind, static=self%static, layered=layered_, &
        allow_static=self%allow_static)], &
        timestamp=timestamp, grid=grid, grid_init_var=grid_init_var, tol=1.0e-5_dp)
    else
      call self%ds%init(path=self%path, vars=[var(name=trim(self%name), kind=kind, static=self%static, layered=layered_, &
        allow_static=self%allow_static)], &
        timestamp=timestamp, grid=grid, tol=1.0e-5_dp)
    end if
    ! TODO: make tol configurable in case of single precision variables with small values
    self%static = self%ds%static
    self%stepping = merge(no_time, self%ds%timestep, self%static)
    self%var_id = self%ds%var_index(self%name)
  end subroutine input_var_open_dataset

  !> \brief Check whether a coupled variable is properly initialized and within time bounds.
  function input_var_check_couple_status(self) result(is_coupled)
    class(input_var_abc), intent(inout) :: self
    logical :: is_coupled !< whether this variable is coupled and should not be updated
    is_coupled = self%coupled
    ! if (is_coupled) then
    !   if (self%static) then
    !     if (.not.allocated(self%cache)) then
    !       log_fatal(*) "Coupled static input variable not initialized: ", trim(self%name)
    !       error stop 1
    !     end if
    !   else if (time <= self%chunk_time_start .or. time > self%chunk_time_end) then
    !     log_fatal(*) "Coupled input variable time out of valid bounds: ", trim(self%name)
    !     error stop 1
    !   end if
    ! end if
  end function input_var_check_couple_status

  !> \brief Set a mask for this input variable from a coupler.
  subroutine input_var_set_mask(self, mask)
    class(input_var_abc), intent(inout) :: self
    logical, intent(in) :: mask(:,:) !< mask to set for this variable
    self%mask = mask
    scope_info(s,*) "Set grid mask from input variable '", trim(self%name), "'"
  end subroutine input_var_set_mask

  !> \brief Update a single double precision input variable.
  subroutine input_var_dp_update(self, exchange_var, time, chunking, end_time)
    class(input_var_dp), intent(inout), target :: self
    real(dp), intent(inout), pointer :: exchange_var(:) !< exchange variable to update
    type(datetime), intent(in) :: time !< current model time
    integer(i4), intent(in) :: chunking !< chunking configuration (0 single read, -1 daily, -2 monthly, -3 yearly, >0 every n hours)
    type(datetime), intent(in) :: end_time !< end time of the simulation (for chunking limits)
    if (.not.self%provided) return
    if (self%check_couple_status()) return ! coupled variables are not updated here
    if (self%static) return ! static variables do not need to be updated
    ! check chunking
    if (chunking /= 1_i4) then
      ! update chunk
      if (time > self%chunk_time_end) then
        call update_time_frame(self%chunk_time_start, self%chunk_time_end, chunking, end_time)
        scope_debug(s,*) "Read new chunk for '", self%name, "': ", self%chunk_time_start%str(), " to ", self%chunk_time_end%str()
        nullify(exchange_var)
        if (allocated(self%cache)) deallocate(self%cache)
        call self%ds%read_chunk(self%var_id, self%cache, self%chunk_time_start, self%chunk_time_end)
        self%offset = self%ds%time_index(self%chunk_time_start)
      end if
      exchange_var => self%cache(:, self%ds%time_index(time) - self%offset)
    else
      if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, 1))
      select case (self%stepping)
        case (1_i4)
          call self%ds%read(self%var_id, self%cache(:, 1), time)
        case (daily, monthly, yearly)
          if (time > self%chunk_time_end) then
            self%chunk_time_start = self%chunk_time_end
            self%chunk_time_end = self%ds%times(self%ds%time_index(time))
            call self%ds%read(self%var_id, self%cache(:, 1), time)
          end if
        case default
          if (self%stepping > 1_i4) then
            if (time > self%chunk_time_end) then
              self%chunk_time_start = self%chunk_time_end
              self%chunk_time_end = self%ds%times(self%ds%time_index(time))
              call self%ds%read(self%var_id, self%cache(:, 1), time)
            end if
          else
            log_fatal(*) "Input: unsupported temporal stepping for variable '", trim(self%name), "': ", n2s(self%stepping)
            error stop 1
          end if
      end select
      exchange_var => self%cache(:, 1)
    end if
  end subroutine input_var_dp_update

  !> \brief Update a single integer input variable.
  subroutine input_var_i4_update(self, exchange_var, time, chunking, end_time)
    class(input_var_i4), intent(inout), target :: self
    integer(i4), intent(inout), pointer :: exchange_var(:) !< exchange variable to update
    type(datetime), intent(in) :: time !< current model time
    integer(i4), intent(in) :: chunking !< chunking configuration (0 single read, -1 daily, -2 monthly, -3 yearly, >0 every n hours)
    type(datetime), intent(in) :: end_time !< end time of the simulation (for chunking limits)
    if (.not.self%provided) return
    if (self%check_couple_status()) return ! coupled variables are not updated here
    if (self%static) return ! static variables do not need to be updated
    ! check chunking
    if (chunking /= 1_i4) then
      ! update chunk
      if (time > self%chunk_time_end) then
        call update_time_frame(self%chunk_time_start, self%chunk_time_end, chunking, end_time)
        scope_debug(s,*) "Read new chunk for '", self%name, "': ", self%chunk_time_start%str(), " to ", self%chunk_time_end%str()
        nullify(exchange_var)
        if (allocated(self%cache)) deallocate(self%cache)
        call self%ds%read_chunk(self%var_id, self%cache, self%chunk_time_start, self%chunk_time_end)
        self%offset = self%ds%time_index(self%chunk_time_start)
      end if
      exchange_var => self%cache(:, self%ds%time_index(time) - self%offset)
    else
      if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, 1))
      select case (self%stepping)
        case (1_i4)
          call self%ds%read(self%var_id, self%cache(:, 1), time)
        case (daily, monthly, yearly)
          if (time > self%chunk_time_end) then
            self%chunk_time_start = self%chunk_time_end
            self%chunk_time_end = self%ds%times(self%ds%time_index(time))
            call self%ds%read(self%var_id, self%cache(:, 1), time)
          end if
        case default
          if (self%stepping > 1_i4) then
            if (time > self%chunk_time_end) then
              self%chunk_time_start = self%chunk_time_end
              self%chunk_time_end = self%ds%times(self%ds%time_index(time))
              call self%ds%read(self%var_id, self%cache(:, 1), time)
            end if
          else
            log_fatal(*) "Input: unsupported temporal stepping for variable '", trim(self%name), "': ", n2s(self%stepping)
            error stop 1
          end if
      end select
      exchange_var => self%cache(:, 1)
    end if
  end subroutine input_var_i4_update

  !> \brief Update a single 16-bit integer input variable.
  subroutine input_var_i2_update(self, exchange_var, time, chunking, end_time)
    class(input_var_i2), intent(inout), target :: self
    integer(i2), intent(inout), pointer :: exchange_var(:) !< exchange variable to update
    type(datetime), intent(in) :: time !< current model time
    integer(i4), intent(in) :: chunking !< chunking configuration (0 single read, -1 daily, -2 monthly, -3 yearly, >0 every n hours)
    type(datetime), intent(in) :: end_time !< end time of the simulation (for chunking limits)
    if (.not.self%provided) return
    if (self%check_couple_status()) return ! coupled variables are not updated here
    if (self%static) return ! static variables do not need to be updated
    ! check chunking
    if (chunking /= 1_i4) then
      ! update chunk
      if (time > self%chunk_time_end) then
        call update_time_frame(self%chunk_time_start, self%chunk_time_end, chunking, end_time)
        scope_debug(s,*) "Read new chunk for '", self%name, "': ", self%chunk_time_start%str(), " to ", self%chunk_time_end%str()
        nullify(exchange_var)
        if (allocated(self%cache)) deallocate(self%cache)
        call self%ds%read_chunk(self%var_id, self%cache, self%chunk_time_start, self%chunk_time_end)
        self%offset = self%ds%time_index(self%chunk_time_start)
      end if
      exchange_var => self%cache(:, self%ds%time_index(time) - self%offset)
    else
      if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, 1))
      select case (self%stepping)
        case (1_i4)
          call self%ds%read(self%var_id, self%cache(:, 1), time)
        case (daily, monthly, yearly)
          if (time > self%chunk_time_end) then
            self%chunk_time_start = self%chunk_time_end
            self%chunk_time_end = self%ds%times(self%ds%time_index(time))
            call self%ds%read(self%var_id, self%cache(:, 1), time)
          end if
        case default
          if (self%stepping > 1_i4) then
            if (time > self%chunk_time_end) then
              self%chunk_time_start = self%chunk_time_end
              self%chunk_time_end = self%ds%times(self%ds%time_index(time))
              call self%ds%read(self%var_id, self%cache(:, 1), time)
            end if
          else
            log_fatal(*) "Input: unsupported temporal stepping for variable '", trim(self%name), "': ", n2s(self%stepping)
            error stop 1
          end if
      end select
      exchange_var => self%cache(:, 1)
    end if
  end subroutine input_var_i2_update

  !> \brief Update a single 2D double precision input variable.
  subroutine input_var2d_dp_update(self, exchange_var, time, chunking, end_time)
    class(input_var2d_dp), intent(inout), target :: self
    real(dp), intent(inout), pointer :: exchange_var(:,:) !< exchange variable to update
    type(datetime), intent(in) :: time !< current model time
    integer(i4), intent(in) :: chunking !< chunking configuration (0 single read, -1 daily, -2 monthly, -3 yearly, >0 every n hours)
    type(datetime), intent(in) :: end_time !< end time of the simulation (for chunking limits)
    if (.not.self%provided) return
    if (self%check_couple_status()) return ! coupled variables are not updated here
    if (self%static) return ! static variables do not need to be updated
    ! check chunking
    if (chunking /= 1_i4) then
      ! update chunk
      if (time > self%chunk_time_end) then
        call update_time_frame(self%chunk_time_start, self%chunk_time_end, chunking, end_time)
        scope_debug(s,*) "Read new chunk for '", self%name, "': ", self%chunk_time_start%str(), " to ", self%chunk_time_end%str()
        nullify(exchange_var)
        if (allocated(self%cache)) deallocate(self%cache)
        call self%ds%read_chunk_layered(self%var_id, self%cache, self%chunk_time_start, self%chunk_time_end)
        self%offset = self%ds%time_index(self%chunk_time_start)
      end if
      exchange_var => self%cache(:, :, self%ds%time_index(time) - self%offset)
    else
      if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, self%ds%nlayers, 1))
      select case (self%stepping)
        case (1_i4)
          call self%ds%read_layered(self%var_id, self%cache(:, :, 1), time)
        case (daily, monthly, yearly)
          if (time > self%chunk_time_end) then
            self%chunk_time_start = self%chunk_time_end
            self%chunk_time_end = self%ds%times(self%ds%time_index(time))
            call self%ds%read_layered(self%var_id, self%cache(:, :, 1), time)
          end if
        case default
          if (self%stepping > 1_i4) then
            if (time > self%chunk_time_end) then
              self%chunk_time_start = self%chunk_time_end
              self%chunk_time_end = self%ds%times(self%ds%time_index(time))
              call self%ds%read_layered(self%var_id, self%cache(:, :, 1), time)
            end if
          else
            log_fatal(*) "Input: unsupported temporal stepping for variable '", trim(self%name), "': ", n2s(self%stepping)
            error stop 1
          end if
      end select
      exchange_var => self%cache(:, :, 1)  ! maybe not needed to re-associate every time
    end if
  end subroutine input_var2d_dp_update

  !> \brief Update a single 2D integer input variable.
  subroutine input_var2d_i4_update(self, exchange_var, time, chunking, end_time)
    class(input_var2d_i4), intent(inout), target :: self
    integer(i4), intent(inout), pointer :: exchange_var(:,:) !< exchange variable to update
    type(datetime), intent(in) :: time !< current model time
    integer(i4), intent(in) :: chunking !< chunking configuration (0 single read, -1 daily, -2 monthly, -3 yearly, >0 every n hours)
    type(datetime), intent(in) :: end_time !< end time of the simulation (for chunking limits)
    if (.not.self%provided) return
    if (self%check_couple_status()) return ! coupled variables are not updated here
    if (self%static) return ! static variables do not need to be updated
    ! check chunking
    if (chunking /= 1_i4) then
      ! update chunk
      if (time > self%chunk_time_end) then
        call update_time_frame(self%chunk_time_start, self%chunk_time_end, chunking, end_time)
        scope_debug(s,*) "Read new chunk for '", self%name, "': ", self%chunk_time_start%str(), " to ", self%chunk_time_end%str()
        nullify(exchange_var)
        if (allocated(self%cache)) deallocate(self%cache)
        call self%ds%read_chunk_layered(self%var_id, self%cache, self%chunk_time_start, self%chunk_time_end)
        self%offset = self%ds%time_index(self%chunk_time_start)
      end if
      exchange_var => self%cache(:, :, self%ds%time_index(time) - self%offset)
    else
      if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, self%ds%nlayers, 1))
      select case (self%stepping)
        case (1_i4)
          call self%ds%read_layered(self%var_id, self%cache(:, :, 1), time)
        case (daily, monthly, yearly)
          if (time > self%chunk_time_end) then
            self%chunk_time_start = self%chunk_time_end
            self%chunk_time_end = self%ds%times(self%ds%time_index(time))
            call self%ds%read_layered(self%var_id, self%cache(:, :, 1), time)
          end if
        case default
          if (self%stepping > 1_i4) then
            if (time > self%chunk_time_end) then
              self%chunk_time_start = self%chunk_time_end
              self%chunk_time_end = self%ds%times(self%ds%time_index(time))
              call self%ds%read_layered(self%var_id, self%cache(:, :, 1), time)
            end if
          else
            log_fatal(*) "Input: unsupported temporal stepping for variable '", trim(self%name), "': ", n2s(self%stepping)
            error stop 1
          end if
      end select
      exchange_var => self%cache(:, :, 1)  ! maybe not needed to re-associate every time
    end if
  end subroutine input_var2d_i4_update

  !> \brief Read a single static double precision input variable and close the dataset.
  subroutine input_var_dp_read_static(self)
    class(input_var_dp), intent(inout), target :: self
    real(dp), allocatable :: data2d(:, :)
    if (.not.self%provided) return
    if (self%coupled) return
    if (.not.self%static) return
    scope_info(s,*) "Read static input variable '", trim(self%name), "'"
    if (self%is_ascii()) then
      if (.not.associated(self%grid)) then
        log_fatal(*) "Input: grid not connected for static ASCII variable: ", trim(self%name)
        error stop 1
      end if
      if (.not.allocated(self%cache)) allocate(self%cache(self%grid%ncells, 1))
      call self%grid%read_data(self%path, data2d)
      call self%grid%pack_into(data2d, self%cache(:, 1))
      if (allocated(data2d)) deallocate(data2d)
      return
    end if
    if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, 1))
    call self%ds%read(self%var_id, self%cache(:, 1))
    call self%ds%close()
  end subroutine input_var_dp_read_static

  !> \brief Read a single static integer input variable and close the dataset.
  subroutine input_var_i4_read_static(self)
    class(input_var_i4), intent(inout), target :: self
    integer(i4), allocatable :: data2d(:, :)
    if (.not.self%provided) return
    if (self%coupled) return
    if (.not.self%static) return
    scope_info(s,*) "Read static input variable '", trim(self%name), "'"
    if (self%is_ascii()) then
      if (.not.associated(self%grid)) then
        log_fatal(*) "Input: grid not connected for static ASCII variable: ", trim(self%name)
        error stop 1
      end if
      if (.not.allocated(self%cache)) allocate(self%cache(self%grid%ncells, 1))
      call self%grid%read_data(self%path, data2d)
      call self%grid%pack_into(data2d, self%cache(:, 1))
      if (allocated(data2d)) deallocate(data2d)
      return
    end if
    if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, 1))
    call self%ds%read(self%var_id, self%cache(:, 1))
    call self%ds%close()
  end subroutine input_var_i4_read_static

  !> \brief Read a single static 16-bit integer input variable and close the dataset.
  subroutine input_var_i2_read_static(self)
    class(input_var_i2), intent(inout), target :: self
    integer(i4), parameter :: i2_min = -int(huge(0_i2), i4) - 1_i4
    integer(i4), parameter :: i2_max = int(huge(0_i2), i4)
    integer(i4), allocatable :: data2d(:, :), data(:)
    if (.not.self%provided) return
    if (self%coupled) return
    if (.not.self%static) return
    scope_info(s,*) "Read static input variable '", trim(self%name), "'"
    if (self%is_ascii()) then
      if (.not.associated(self%grid)) then
        log_fatal(*) "Input: grid not connected for static ASCII variable: ", trim(self%name)
        error stop 1
      end if
      if (.not.allocated(self%cache)) allocate(self%cache(self%grid%ncells, 1))
      call self%grid%read_data(self%path, data2d)
      allocate(data(self%grid%ncells))
      call self%grid%pack_into(data2d, data)
      if (any(data < i2_min .or. data > i2_max)) then
        log_fatal(*) "Input: ASCII data outside the i2 range for variable: ", trim(self%name)
        error stop 1
      end if
      self%cache(:, 1) = int(data, i2)
      if (allocated(data)) deallocate(data)
      if (allocated(data2d)) deallocate(data2d)
      return
    end if
    if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, 1))
    call self%ds%read(self%var_id, self%cache(:, 1))
    call self%ds%close()
  end subroutine input_var_i2_read_static

  !> \brief Read a single static 2D double precision input variable and close the dataset.
  subroutine input_var2d_dp_read_static(self)
    class(input_var2d_dp), intent(inout), target :: self
    if (.not.self%provided) return
    if (self%coupled) return
    if (.not.self%static) return
    scope_info(s,*) "Read static input variable '", trim(self%name), "'"
    if (self%is_ascii()) then
      log_fatal(*) "Input: layered ASCII input is not supported for variable: ", trim(self%name)
      error stop 1
    end if
    if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, self%ds%nlayers, 1))
    call self%ds%read_layered(self%var_id, self%cache(:, :, 1))
    call self%ds%close()
  end subroutine input_var2d_dp_read_static

  !> \brief Read a single static 2D integer input variable and close the dataset.
  subroutine input_var2d_i4_read_static(self)
    class(input_var2d_i4), intent(inout), target :: self
    if (.not.self%provided) return
    if (self%coupled) return
    if (.not.self%static) return
    scope_info(s,*) "Read static input variable '", trim(self%name), "'"
    if (self%is_ascii()) then
      log_fatal(*) "Input: layered ASCII input is not supported for variable: ", trim(self%name)
      error stop 1
    end if
    if (.not.allocated(self%cache)) allocate(self%cache(self%ds%grid%ncells, self%ds%nlayers, 1))
    call self%ds%read_layered(self%var_id, self%cache(:, :, 1))
    call self%ds%close()
  end subroutine input_var2d_i4_read_static

  !> \brief Reset the time frame for chunked reading.
  subroutine input_var_dp_reset_time(self, chunking, start_time)
    class(input_var_dp), intent(inout) :: self
    integer(i4), intent(in) :: chunking !< chunking mode
    type(datetime), intent(in) :: start_time !< start time of the simulation
    if (.not.self%provided) return
    if (self%static) return
    ! if read once and cache already allocated, do nothing
    if (.not.self%coupled .and. chunking == 0_i4 .and. allocated(self%cache)) return
    self%chunk_time_start = start_time
    self%chunk_time_end = start_time ! trigger read on first update or require coupler set
  end subroutine input_var_dp_reset_time

  !> \brief Reset the time frame for chunked reading.
  subroutine input_var_i4_reset_time(self, chunking, start_time)
    class(input_var_i4), intent(inout) :: self
    integer(i4), intent(in) :: chunking !< chunking mode
    type(datetime), intent(in) :: start_time !< start time of the simulation
    if (.not.self%provided) return
    if (self%static) return
    ! if read once and cache already allocated, do nothing
    if (.not.self%coupled .and. chunking == 0_i4 .and. allocated(self%cache)) return
    self%chunk_time_start = start_time
    self%chunk_time_end = start_time ! trigger read on first update or require coupler set
  end subroutine input_var_i4_reset_time

  !> \brief Reset the time frame for chunked reading.
  subroutine input_var_i2_reset_time(self, chunking, start_time)
    class(input_var_i2), intent(inout) :: self
    integer(i4), intent(in) :: chunking !< chunking mode
    type(datetime), intent(in) :: start_time !< start time of the simulation
    if (.not.self%provided) return
    if (self%static) return
    ! if read once and cache already allocated, do nothing
    if (.not.self%coupled .and. chunking == 0_i4 .and. allocated(self%cache)) return
    self%chunk_time_start = start_time
    self%chunk_time_end = start_time ! trigger read on first update or require coupler set
  end subroutine input_var_i2_reset_time

  !> \brief Reset the time frame for chunked reading.
  subroutine input_var2d_dp_reset_time(self, chunking, start_time)
    class(input_var2d_dp), intent(inout) :: self
    integer(i4), intent(in) :: chunking !< chunking mode
    type(datetime), intent(in) :: start_time !< start time of the simulation
    if (.not.self%provided) return
    if (self%static) return
    ! if read once and cache already allocated, do nothing
    if (.not.self%coupled .and. chunking == 0_i4 .and. allocated(self%cache)) return
    self%chunk_time_start = start_time
    self%chunk_time_end = start_time ! trigger read on first update or require coupler set
  end subroutine input_var2d_dp_reset_time

  !> \brief Reset the time frame for chunked reading.
  subroutine input_var2d_i4_reset_time(self, chunking, start_time)
    class(input_var2d_i4), intent(inout) :: self
    integer(i4), intent(in) :: chunking !< chunking mode
    type(datetime), intent(in) :: start_time !< start time of the simulation
    if (.not.self%provided) return
    if (self%static) return
    ! if read once and cache already allocated, do nothing
    if (.not.self%coupled .and. chunking == 0_i4 .and. allocated(self%cache)) return
    self%chunk_time_start = start_time
    self%chunk_time_end = start_time ! trigger read on first update or require coupler set
  end subroutine input_var2d_i4_reset_time

  !> \brief Initialize the input configuration.
  subroutine input_config_read(self, file)
    use nml_helper, only: NML_OK
    ! input/output variables
    class(input_config_t), target, intent(inout) :: self
    character(*), intent(in) :: file !< file containing the namelists
    character(1024) :: errmsg
    integer :: status
    log_info(*) "Read input configuration from file: ", trim(file)
    status = self%input%from_file(file=file, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Input: Error reading input configuration from file: ", trim(file), ", with error: ", trim(errmsg)
      error stop 1
    end if
    ! TODO: read coupling configuration
  end subroutine input_config_read

  !> \brief Set runtime dimensions for generated input namelists.
  subroutine input_set_dims(self)
    class(input_t), target, intent(inout) :: self
    character(1024) :: errmsg
    integer :: status

    status = self%config%input%set_dims(n_domains=self%exchange%nml_n_domains, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Input: Error setting input dimensions: ", trim(errmsg)
      error stop 1
    end if
    status = self%config%coupling%set_dims(n_domains=self%exchange%nml_n_domains, errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Input: Error setting coupling dimensions: ", trim(errmsg)
      error stop 1
    end if
  end subroutine input_set_dims

  !> \brief Configure the Input container.
  !> \details Read configuration from file and initialize input variables.
  !! If file is not given, the configuration is assumed to be already set from a coupler,
  !! which is caught by the "is_configured" check.
  subroutine input_configure(self, file)
    class(input_t), target, intent(inout) :: self
    character(*), intent(in), optional :: file !< file containing the namelists
    integer :: status
    character(1024) :: errmsg
    character(:), allocatable :: path
    integer(i4) :: id(1)
    log_info(*) "Configure Input"
    if (present(file)) then
      path = self%exchange%get_path(file)
      call self%config%read(path)
    end if
    if (.not.self%config%input%is_configured) then
      log_fatal(*) "Input configuration not set."
      error stop 1
    end if
    status = self%config%input%is_valid(errmsg=errmsg)
    if (status /= NML_OK) then
      log_fatal(*) "Input configuration not valid: ", trim(errmsg)
      error stop 1
    end if

    ! initialize general inputs settings
    id(1) = self%exchange%nml_domain_id ! domain ID to read correct namelist entries
    self%chunking = self%config%input%chunking(id(1))
    self%time_stamp_location = self%config%input%time_stamp_location(id(1))
    self%morph_latlon = self%config%input%morph_latlon(id(1))

    ! morph mask (morph_mask_var by default "mask")
    status = self%config%input%is_set("morph_mask_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%morph_mask%init( &
        path=self%exchange%get_path(self%config%input%morph_mask_path(id(1))), name=self%config%input%morph_mask_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
    end if

    ! DEM (dem_var by default "dem")
    status = self%config%input%is_set("dem_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%dem%init( &
        path=self%exchange%get_path(self%config%input%dem_path(id(1))), name=self%config%input%dem_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
      self%exchange%dem%provided = .true. ! mark as provided in exchange
    end if

    ! slope (slope_var by default "slope")
    status = self%config%input%is_set("slope_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%slope%init( &
        path=self%exchange%get_path(self%config%input%slope_path(id(1))), name=self%config%input%slope_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
      self%exchange%slope%provided = .true. ! mark as provided in exchange
    end if

    ! aspect (aspect_var by default "aspect")
    status = self%config%input%is_set("aspect_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%aspect%init( &
        path=self%exchange%get_path(self%config%input%aspect_path(id(1))), name=self%config%input%aspect_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
      self%exchange%aspect%provided = .true. ! mark as provided in exchange
    end if

    ! flow direction (fdir_var by default "fdir")
    status = self%config%input%is_set("fdir_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%fdir%init( &
        path=self%exchange%get_path(self%config%input%fdir_path(id(1))), name=self%config%input%fdir_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
    end if

    ! flow accumulation (facc_var by default "facc")
    status = self%config%input%is_set("facc_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%facc%init( &
        path=self%exchange%get_path(self%config%input%facc_path(id(1))), name=self%config%input%facc_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
      self%exchange%facc%provided = .true. ! mark as provided in exchange
    end if

    ! geology class (geo_class_var by default "geology_class")
    status = self%config%input%is_set("geo_class_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%geo_class%init( &
        path=self%exchange%get_path(self%config%input%geo_class_path(id(1))), name=self%config%input%geo_class_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
      self%exchange%geo_unit%provided = .true. ! mark as provided in exchange
    end if

    ! soil class (soil_class_var by default "soil_class")
    status = self%config%input%is_set("soil_class_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%soil_class%init( &
        path=self%exchange%get_path(self%config%input%soil_class_path(id(1))), name=self%config%input%soil_class_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
      self%exchange%soil_id%provided = .true. ! mark as provided in exchange
    end if

    ! soil horizon class (uses soil_class_var by default)
    status = self%config%input%is_set("soil_horizon_class_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%soil_horizon_class%init( &
        path=self%exchange%get_path(self%config%input%soil_horizon_class_path(id(1))), name=self%config%input%soil_class_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
      self%exchange%soil_id%provided = .true. ! mark as provided in exchange
    end if

    if (self%soil_class%provided .and. self%soil_horizon_class%provided) then
      log_warn(*) "Input: both soil_class_path and soil_horizon_class_path are configured. Prefer layered soil horizon class."
    end if

    ! LAI class (lai_class_var by default "LAI_class")
    status = self%config%input%is_set("lai_class_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%lai_class%init( &
        path=self%exchange%get_path(self%config%input%lai_class_path(id(1))), name=self%config%input%lai_class_var(id(1)), &
        static=.true., morph_latlon=self%morph_latlon)
      self%exchange%lai_class%provided = .true. ! mark as provided in exchange
    end if

    ! meteorological mask (meteo_mask_var by default "mask")
    status = self%config%input%is_set("meteo_mask_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%meteo_mask%init( &
        path=self%exchange%get_path(self%config%input%meteo_mask_path(id(1))), name=self%config%input%meteo_mask_var(id(1)), allow_static=.true.)
    end if

    ! raw precipitation on level2
    status = self%config%input%is_set("pre_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%pre%init(path=self%exchange%get_path(self%config%input%pre_path(id(1))), name=self%config%input%pre_var(id(1)))
      self%exchange%raw_pre%provided = .true.
    end if

    ! raw PET on level2
    status = self%config%input%is_set("pet_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%pet%init(path=self%exchange%get_path(self%config%input%pet_path(id(1))), name=self%config%input%pet_var(id(1)))
      self%exchange%raw_pet%provided = .true.
    end if

    ! raw temperature on level2
    status = self%config%input%is_set("temp_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%temp%init(path=self%exchange%get_path(self%config%input%temp_path(id(1))), name=self%config%input%temp_var(id(1)))
      self%exchange%raw_temp%provided = .true.
    end if

    ! raw annual mean temperature on level2
    status = self%config%input%is_set("tann_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%tann%init(path=self%exchange%get_path(self%config%input%tann_path(id(1))), name=self%config%input%tann_var(id(1)), allow_static=.true.)
      self%exchange%raw_tann%provided = .true.
    end if

    ! raw minimum temperature on level2
    status = self%config%input%is_set("tmin_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%tmin%init(path=self%exchange%get_path(self%config%input%tmin_path(id(1))), name=self%config%input%tmin_var(id(1)))
      self%exchange%raw_tmin%provided = .true.
    end if

    ! raw maximum temperature on level2
    status = self%config%input%is_set("tmax_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%tmax%init(path=self%exchange%get_path(self%config%input%tmax_path(id(1))), name=self%config%input%tmax_var(id(1)))
      self%exchange%raw_tmax%provided = .true.
    end if

    ! raw short-wave radiation on level2
    status = self%config%input%is_set("ssrd_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%ssrd%init(path=self%exchange%get_path(self%config%input%ssrd_path(id(1))), name=self%config%input%ssrd_var(id(1)))
      self%exchange%raw_ssrd%provided = .true.
    end if

    ! raw long-wave radiation on level2
    status = self%config%input%is_set("strd_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%strd%init(path=self%exchange%get_path(self%config%input%strd_path(id(1))), name=self%config%input%strd_var(id(1)))
      self%exchange%raw_strd%provided = .true.
    end if

    ! raw net radiation on level2
    status = self%config%input%is_set("netrad_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%netrad%init(path=self%exchange%get_path(self%config%input%netrad_path(id(1))), name=self%config%input%netrad_var(id(1)))
      self%exchange%raw_netrad%provided = .true.
    end if

    ! raw vapor pressure on level2
    status = self%config%input%is_set("eabs_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%eabs%init(path=self%exchange%get_path(self%config%input%eabs_path(id(1))), name=self%config%input%eabs_var(id(1)))
      self%exchange%raw_eabs%provided = .true.
    end if

    ! raw wind speed on level2
    status = self%config%input%is_set("wind_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%wind%init(path=self%exchange%get_path(self%config%input%wind_path(id(1))), name=self%config%input%wind_var(id(1)))
      self%exchange%raw_wind%provided = .true.
    end if

    ! hydro mask (hydro_mask_var by default "mask")
    status = self%config%input%is_set("hydro_mask_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%hydro_mask%init( &
        path=self%exchange%get_path(self%config%input%hydro_mask_path(id(1))), name=self%config%input%hydro_mask_var(id(1)), static=.true.)
    end if

    ! runoff (runoff_var by default "runoff")
    status = self%config%input%is_set("runoff_path", idx=id, errmsg=errmsg)
    if (status == NML_OK) then
      call self%runoff%init( &
        path=self%exchange%get_path(self%config%input%runoff_path(id(1))), name=self%config%input%runoff_var(id(1)))
      self%exchange%runoff_total%provided = .true. ! mark as provided in exchange
    end if
  end subroutine input_configure

  !> \brief Connect the Input container.
  !> \details Open datasets and initialize grids if needed.
  !! Principle: first come first serve for grid initialization.
  !! All static variables are read here and available afterwards for other containers.
  subroutine input_connect(self)
    class(input_t), target, intent(inout) :: self
    logical :: init_grid
    integer(i4) :: ts
    log_info(*) "Connect Input"
    ts = self%time_stamp_location

    ! Flow direction establishes the full level-0 grid and river before all other morphology inputs.
    if (self%fdir%coupled) then
      log_fatal(*) "Input: coupled flow direction and level-0 grid initialization is not implemented."
      error stop 1
    else if (self%fdir%provided) then
      call self%read_fdir_file()
    end if

    ! morph mask
    if (self%morph_mask%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: morph mask is coupled... not yet implemented"
      stop 1
    else if (self%morph_mask%provided) then
      init_grid = need_grid(self%tgt_level0, self%exchange%level0) ! associate grid if not yet done
      call self%morph_mask%open_dataset(kind="i4", timestamp=ts, grid=self%exchange%level0, init_grid=init_grid)
      if (.not.self%morph_mask%is_ascii()) call self%morph_mask%ds%close() ! morph mask not read, since we only need the grid information
      call self%morph_mask%set_mask(self%exchange%level0%mask) ! only done for masks
    end if

    ! DEM
    if (self%dem%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: DEM is coupled... not yet implemented"
      stop 1
    else if (self%dem%provided) then
      init_grid = need_grid(self%tgt_level0, self%exchange%level0) ! associate grid if not yet done
      call self%dem%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level0, init_grid=init_grid)
      call self%dem%read_static()
      if (self%owns_river_l0) call self%river_l0%set_elevation(self%dem%cache(:, 1))
    end if

    ! Lake delineation needs the fdir-derived river and the full packed DEM.
    if (self%config%input%is_set("lake_definition_path", idx=[self%exchange%nml_domain_id]) == NML_OK) &
      call self%read_lake_specification()
    if (.not.associated(self%exchange%level0_land) .and. associated(self%exchange%level0)) &
      self%exchange%level0_land => self%exchange%level0
    if (self%dem%provided) then
      if (self%owns_river_l0 .and. .not.associated(self%exchange%level0_lake)) then
        ! The full DEM belongs to the fdir-derived river; no land-packed view is needed without lakes.
        self%exchange%dem%data => self%river_l0%node_elevation
        deallocate(self%dem%cache)
      else
        self%exchange%dem%data => self%dem%cache(:, 1)
      end if
    end if

    ! slope
    if (self%slope%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: slope is coupled... not yet implemented"
      stop 1
    else if (self%slope%provided) then
      if (self%owns_river_l0) then
        init_grid = need_grid(self%tgt_level0, self%exchange%level0)
        call self%slope%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level0, init_grid=init_grid)
      else
        init_grid = need_level0_land_grid(self%tgt_level0, self%exchange%level0, self%exchange%level0_land)
        call self%slope%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level0_land, init_grid=init_grid)
      end if
      call self%slope%read_static()
      if (self%owns_river_l0) then
        call self%river_l0%set_link_slope(self%slope%cache(:, 1))
        if (associated(self%exchange%level0_lake)) &
          call repack_l0_dp_cache(self%exchange%level0, self%exchange%level0_land, self%slope%cache, "slope")
        self%slope%grid => self%exchange%level0_land
      end if
      self%exchange%slope%data => self%slope%cache(:, 1) ! associate exchange variable to input cache
    end if

    ! aspect
    if (self%aspect%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: aspect is coupled... not yet implemented"
      stop 1
    else if (self%aspect%provided) then
      init_grid = need_level0_land_grid(self%tgt_level0, self%exchange%level0, self%exchange%level0_land)
      call self%aspect%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level0_land, init_grid=init_grid)
      call self%aspect%read_static()
      self%exchange%aspect%data => self%aspect%cache(:, 1) ! associate exchange variable to input cache
    end if

    ! flow accumulation
    if (self%facc%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: flow accumulation is coupled... not yet implemented"
      stop 1
    else if (self%facc%provided) then
      init_grid = need_grid(self%tgt_level0, self%exchange%level0) ! associate grid if not yet done
      call self%facc%open_dataset(kind="i4", timestamp=ts, grid=self%exchange%level0, init_grid=init_grid)
      call self%facc%read_static()
      self%exchange%facc%data => self%facc%cache(:, 1) ! associate exchange variable to input cache
    end if

    ! geology class
    if (self%geo_class%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: geology class is coupled... not yet implemented"
      stop 1
    else if (self%geo_class%provided) then
      init_grid = need_level0_land_grid(self%tgt_level0, self%exchange%level0, self%exchange%level0_land)
      call self%geo_class%open_dataset(kind="i4", timestamp=ts, grid=self%exchange%level0_land, init_grid=init_grid)
      call self%geo_class%read_static()
      self%exchange%geo_unit%data => self%geo_class%cache(:, 1) ! associate exchange variable to input cache
    end if

    ! soil class
    if (self%soil_class%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: soil class is coupled... not yet implemented"
      stop 1
    else if (self%soil_class%provided) then
      init_grid = need_level0_land_grid(self%tgt_level0, self%exchange%level0, self%exchange%level0_land)
      call self%soil_class%open_dataset(kind="i4", timestamp=ts, grid=self%exchange%level0_land, init_grid=init_grid)
      call self%soil_class%read_static()
      if (.not.self%soil_horizon_class%provided) then
        if (allocated(self%soil_class_one_layer)) deallocate(self%soil_class_one_layer)
        allocate(self%soil_class_one_layer(size(self%soil_class%cache, 1), 1))
        self%soil_class_one_layer(:, 1) = self%soil_class%cache(:, 1)
        self%exchange%soil_id%data => self%soil_class_one_layer ! map single-layer soil class to 2D exchange shape
      end if
    end if

    ! soil horizon class
    if (self%soil_horizon_class%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: soil horizon class is coupled... not yet implemented"
      stop 1
    else if (self%soil_horizon_class%provided) then
      init_grid = need_level0_land_grid(self%tgt_level0, self%exchange%level0, self%exchange%level0_land)
      call self%soil_horizon_class%open_dataset( &
        kind="i4", timestamp=ts, grid=self%exchange%level0_land, init_grid=init_grid, layered=.true.)
      call self%soil_horizon_class%read_static()
      self%exchange%soil_id%data => self%soil_horizon_class%cache(:, :, 1) ! associate exchange variable to input cache
    end if

    ! LAI class
    if (self%lai_class%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: LAI class is coupled... not yet implemented"
      stop 1
    else if (self%lai_class%provided) then
      init_grid = need_level0_land_grid(self%tgt_level0, self%exchange%level0, self%exchange%level0_land)
      call self%lai_class%open_dataset(kind="i4", timestamp=ts, grid=self%exchange%level0_land, init_grid=init_grid)
      call self%lai_class%read_static()
      self%exchange%lai_class%data => self%lai_class%cache(:, 1) ! associate exchange variable to input cache
    end if

    ! meteorological mask
    if (self%meteo_mask%coupled) then
      log_error(*) "Input: meteorological mask is coupled... not yet implemented"
      stop 1
    else if (self%meteo_mask%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%meteo_mask%open_dataset(kind="i4", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      if (.not.self%meteo_mask%is_ascii()) call self%meteo_mask%ds%close()
      call self%meteo_mask%set_mask(self%exchange%level2%mask)
    end if

    ! raw precipitation
    if (self%pre%coupled) then
      log_error(*) "Input: precipitation is coupled... not yet implemented"
      stop 1
    else if (self%pre%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%pre%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%pre, self%exchange%raw_pre)
      if (self%pre%static) then
        call self%pre%read_static()
        self%exchange%raw_pre%data => self%pre%cache(:, 1)
      end if
    end if

    ! raw PET
    if (self%pet%coupled) then
      log_error(*) "Input: PET is coupled... not yet implemented"
      stop 1
    else if (self%pet%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%pet%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%pet, self%exchange%raw_pet)
      if (self%pet%static) then
        call self%pet%read_static()
        self%exchange%raw_pet%data => self%pet%cache(:, 1)
      end if
    end if

    ! raw temperature
    if (self%temp%coupled) then
      log_error(*) "Input: temperature is coupled... not yet implemented"
      stop 1
    else if (self%temp%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%temp%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%temp, self%exchange%raw_temp)
      if (self%temp%static) then
        call self%temp%read_static()
        self%exchange%raw_temp%data => self%temp%cache(:, 1)
      end if
    end if

    ! raw annual mean temperature
    if (self%tann%coupled) then
      log_error(*) "Input: annual mean temperature is coupled... not yet implemented"
      stop 1
    else if (self%tann%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%tann%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%tann, self%exchange%raw_tann)
      if (self%tann%static) then
        call self%tann%read_static()
        self%exchange%raw_tann%data => self%tann%cache(:, 1)
      end if
    end if

    ! raw minimum temperature
    if (self%tmin%coupled) then
      log_error(*) "Input: minimum temperature is coupled... not yet implemented"
      stop 1
    else if (self%tmin%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%tmin%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%tmin, self%exchange%raw_tmin)
      if (self%tmin%static) then
        call self%tmin%read_static()
        self%exchange%raw_tmin%data => self%tmin%cache(:, 1)
      end if
    end if

    ! raw maximum temperature
    if (self%tmax%coupled) then
      log_error(*) "Input: maximum temperature is coupled... not yet implemented"
      stop 1
    else if (self%tmax%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%tmax%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%tmax, self%exchange%raw_tmax)
      if (self%tmax%static) then
        call self%tmax%read_static()
        self%exchange%raw_tmax%data => self%tmax%cache(:, 1)
      end if
    end if

    ! raw short-wave radiation
    if (self%ssrd%coupled) then
      log_error(*) "Input: short-wave radiation is coupled... not yet implemented"
      stop 1
    else if (self%ssrd%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%ssrd%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%ssrd, self%exchange%raw_ssrd)
      if (self%ssrd%static) then
        call self%ssrd%read_static()
        self%exchange%raw_ssrd%data => self%ssrd%cache(:, 1)
      end if
    end if

    ! raw long-wave radiation
    if (self%strd%coupled) then
      log_error(*) "Input: long-wave radiation is coupled... not yet implemented"
      stop 1
    else if (self%strd%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%strd%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%strd, self%exchange%raw_strd)
      if (self%strd%static) then
        call self%strd%read_static()
        self%exchange%raw_strd%data => self%strd%cache(:, 1)
      end if
    end if

    ! raw net radiation
    if (self%netrad%coupled) then
      log_error(*) "Input: net radiation is coupled... not yet implemented"
      stop 1
    else if (self%netrad%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%netrad%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%netrad, self%exchange%raw_netrad)
      if (self%netrad%static) then
        call self%netrad%read_static()
        self%exchange%raw_netrad%data => self%netrad%cache(:, 1)
      end if
    end if

    ! raw vapor pressure
    if (self%eabs%coupled) then
      log_error(*) "Input: vapor pressure is coupled... not yet implemented"
      stop 1
    else if (self%eabs%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%eabs%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%eabs, self%exchange%raw_eabs)
      if (self%eabs%static) then
        call self%eabs%read_static()
        self%exchange%raw_eabs%data => self%eabs%cache(:, 1)
      end if
    end if

    ! raw wind speed
    if (self%wind%coupled) then
      log_error(*) "Input: wind speed is coupled... not yet implemented"
      stop 1
    else if (self%wind%provided) then
      init_grid = need_grid(self%tgt_level2, self%exchange%level2)
      call self%wind%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level2, init_grid=init_grid)
      call sync_input_var_meta(self%wind, self%exchange%raw_wind)
      if (self%wind%static) then
        call self%wind%read_static()
        self%exchange%raw_wind%data => self%wind%cache(:, 1)
      end if
    end if

    ! hydro mask
    if (self%hydro_mask%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: hydro mask is coupled... not yet implemented"
      stop 1
    else if (self%hydro_mask%provided) then
      init_grid = need_grid(self%tgt_level1, self%exchange%level1) ! associate grid if not yet done
      call self%hydro_mask%open_dataset(kind="i4", timestamp=ts, grid=self%exchange%level1, init_grid=init_grid)
      if (.not.self%hydro_mask%is_ascii()) call self%hydro_mask%ds%close() ! hydro mask not read, since we only need the grid information
      call self%hydro_mask%set_mask(self%exchange%level1%mask) ! only done for masks
    end if

    ! runoff
    if (self%runoff%coupled) then
      ! TODO: init grid from coupling namelist if needed
      log_error(*) "Input: runoff is coupled... not yet implemented"
      stop 1
    else if (self%runoff%provided) then
      init_grid = need_grid(self%tgt_level1, self%exchange%level1) ! associate grid if not yet done
      call self%runoff%open_dataset(kind="dp", timestamp=ts, grid=self%exchange%level1, init_grid=init_grid)
      call sync_input_var_meta(self%runoff, self%exchange%runoff_total)
    end if

    if (.not.associated(self%exchange%level0_land) .and. associated(self%exchange%level0)) &
      self%exchange%level0_land => self%exchange%level0

    if (associated(self%exchange%level2)) self%exchange%level2_resolution = self%exchange%level2%cellsize
  end subroutine input_connect

  !> \brief Read file-based flow direction once while constructing the level-0 grid.
  subroutine input_read_fdir_file(self)
    class(input_t), target, intent(inout) :: self
    type(data_t) :: data
    integer(i2), allocatable :: fdir_packed(:)
    logical :: level0_supplied

    level0_supplied = associated(self%exchange%level0)
    data%dtype = "i16"

    scope_info(s,*) "Read flow direction and initialize level-0 grid from file: ", trim(self%fdir%path)
    if (self%fdir%is_ascii()) then
      call self%tgt_level0%from_ascii_file( &
        self%fdir%path, coordsys=merge(spherical, cartesian, self%morph_latlon), data=data)
    else
      call self%tgt_level0%from_netcdf(self%fdir%path, self%fdir%name, tol=1.0e-5_dp, data=data)
    end if
    if (level0_supplied) then
      if (.not.self%exchange%level0%is_matching(self%tgt_level0, tol=1.0e-5_dp)) then
        log_fatal(*) "Input: fdir-derived level-0 grid does not match the supplied level-0 grid."
        error stop 1
      end if
    else
      self%exchange%level0 => self%tgt_level0
    end if

    allocate(fdir_packed(self%exchange%level0%ncells))
    call self%exchange%level0%pack_into(data%data_i2, fdir_packed)
    call data%deallocate()
    call self%build_river_l0(self%exchange%level0, fdir_packed)
    deallocate(fdir_packed)
  end subroutine input_read_fdir_file

  !> \brief Construct and publish the level-0 river from an established grid and packed flow direction.
  !> \details File readers and future couplers should converge on this representation.
  subroutine input_build_river_l0(self, grid, fdir)
    class(input_t), target, intent(inout) :: self
    type(grid_t), pointer, intent(in) :: grid
    integer(i2), intent(in) :: fdir(:)

    if (self%owns_river_l0 .or. associated(self%exchange%river_l0)) then
      log_fatal(*) "Input: level-0 river already constructed before flow-direction handoff."
      error stop 1
    end if
    call self%river_l0%from_fdir(fdir, grid)
    call self%river_l0%calc_order(root=.true.)
    call self%river_l0%calc_facc()
    self%exchange%river_l0 => self%river_l0
    call self%exchange%fdir%publish_local("Input", self%river_l0%fdir, no_time)
    self%owns_river_l0 = .true.
  end subroutine input_build_river_l0

  !> \brief Read static lake metadata, snap outlets to L0 cells, and delineate lake footprints.
  subroutine input_read_lake_specification(self)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    class(input_t), target, intent(inout) :: self
    type(points_input_dataset) :: input
    type(var), allocatable :: vars(:)
    real(dp), allocatable :: coords(:,:)
    integer(i8), allocatable :: outlet_nodes(:)
    character(:), allocatable :: path, name
    integer(i8) :: i, node
    integer(i4) :: ix, iy
    logical :: use_aux

    if (.not.self%owns_river_l0 .or. .not.associated(self%exchange%river_l0)) then
      log_fatal(*) "Input: lake definitions require a file-based level-0 river."
      error stop 1
    end if
    if (.not.self%dem%provided .or. .not.allocated(self%dem%cache)) then
      log_fatal(*) "Input: lake definitions require a static DEM on level 0."
      error stop 1
    end if

    path = self%exchange%get_path(self%config%input%lake_definition_path(self%exchange%nml_domain_id))
    name = trim(self%config%input%lake_level_var(self%exchange%nml_domain_id))
    scope_info(s,*) "Read lake specification from file: ", path
    vars = [var(name=name, long_name="maximum lake level", units="m", static=.true.)]
    call input%init(path, vars=vars, points=self%lake_outlets, points_init_var=name)
    call input%get_ids(self%lake_ids)
    allocate(self%lake_max_levels(self%lake_outlets%n_points))
    call input%read(name, self%lake_max_levels)
    call input%close()

    if (self%lake_outlets%n_points < 1_i8) then
      log_fatal(*) "Input: lake specification contains no lake outlets."
      error stop 1
    end if
    if (size(self%lake_ids, kind=i8) /= self%lake_outlets%n_points .or. &
        size(self%lake_max_levels, kind=i8) /= self%lake_outlets%n_points) then
      log_fatal(*) "Input: lake specification arrays have inconsistent sizes."
      error stop 1
    end if
    if (.not.all(ieee_is_finite(self%lake_outlets%x)) .or. .not.all(ieee_is_finite(self%lake_outlets%y))) then
      log_fatal(*) "Input: lake outlet coordinates must be finite."
      error stop 1
    end if

    coords = self%lake_outlets%coords()
    use_aux = self%lake_outlets%coordsys /= self%exchange%level0%coordsys
    if (use_aux) then
      if (self%lake_outlets%coordsys /= spherical .or. .not.self%exchange%level0%has_aux_vertices()) then
        log_fatal(*) "Input: lake outlet coordinates are incompatible with the level-0 grid."
        error stop 1
      end if
    end if
    outlet_nodes = self%exchange%level0%closest_cell_id(coords, use_aux=use_aux)
    do i = 1_i8, size(outlet_nodes, kind=i8)
      node = outlet_nodes(i)
      if (node < 1_i8 .or. node > self%exchange%level0%ncells) then
        log_fatal(*) "Input: lake outlet cannot be mapped to an active level-0 cell."
        error stop 1
      end if
      ix = self%exchange%level0%cell_ij(node, 1)
      iy = self%exchange%level0%cell_ij(node, 2)
      if (.not.self%exchange%level0%in_cell(ix, iy, coords(i, 1), coords(i, 2), aux=use_aux)) then
        log_fatal(*) "Input: lake outlet lies outside its mapped active level-0 cell."
        error stop 1
      end if
    end do

    call self%river_l0%label_lakes(outlet_nodes, self%lake_ids, self%lake_max_levels)
    call self%build_lake_grids()
    call repack_l0_dp_cache(self%exchange%level0, self%exchange%level0_land, self%dem%cache, "DEM")
    self%dem%grid => self%exchange%level0_land
    self%exchange%lake_points => self%lake_outlets
    call self%exchange%lake_ids%publish_local("Input", self%lake_ids, no_time)
    call self%exchange%lake_max_levels%publish_local("Input", self%lake_max_levels, no_time)
  end subroutine input_read_lake_specification

  !> \brief Construct complementary land and lake grids on the full level-0 geometry.
  subroutine input_build_lake_grids(self)
    class(input_t), target, intent(inout) :: self
    logical, allocatable :: land_packed(:), lake_packed(:)
    logical, allocatable :: new_mask(:,:)

    if (.not.associated(self%exchange%level0)) then
      log_fatal(*) "Input: cannot derive land and lake grids without the full level-0 grid."
      error stop 1
    end if
    call self%river_l0%lake_masks(land_packed, lake_packed)
    allocate(new_mask(self%exchange%level0%nx, self%exchange%level0%ny))
    ! land grid
    call self%exchange%level0%unpack_into(land_packed, new_mask)
    call self%exchange%level0%copy_to(self%tgt_level0_land, mask=new_mask, keep_sub_cell_area=.true.)
    ! lake grid
    call self%exchange%level0%unpack_into(lake_packed, new_mask)
    call self%exchange%level0%copy_to(self%tgt_level0_lake, mask=new_mask, keep_sub_cell_area=.true.)
    deallocate(new_mask)
    if (self%tgt_level0_land%ncells < 1_i8) then
      log_fatal(*) "Input: lake delineation leaves no active land cells."
      error stop 1
    end if
    if (self%tgt_level0_lake%ncells < 1_i8) then
      log_fatal(*) "Input: lake delineation produced no active lake cells."
      error stop 1
    end if
    self%exchange%level0_land => self%tgt_level0_land
    self%exchange%level0_lake => self%tgt_level0_lake
  end subroutine input_build_lake_grids

  !> \brief Initialize the Input container for the model run.
  !> \details Prepare chunked reading if needed.
  subroutine input_initialize(self)
    class(input_t), target, intent(inout) :: self
    log_info(*) "Initialize Input"
    scope_info(s,*) "Initialize input time windows with chunking: ", n2s(self%chunking)
    ! warn about provided but not required variables, since this likely indicates a configuration issue
    if (self%facc%provided .and. .not.self%exchange%facc%required) then
      log_warn(*) "Input: flow accumulation provided but not required. Check your configuration."
    end if
    if (self%dem%provided .and. .not.self%exchange%dem%required .and. &
        .not.associated(self%exchange%level0_lake)) then
      log_warn(*) "Input: DEM provided but not required. Check your configuration."
    end if
    if (self%slope%provided .and. .not.self%exchange%slope%required .and. &
        .not.(self%owns_river_l0 .and. allocated(self%river_l0%link_slope))) then
      log_warn(*) "Input: slope provided but not required. Check your configuration."
    end if
    if (self%aspect%provided .and. .not.self%exchange%aspect%required) then
      log_warn(*) "Input: aspect provided but not required. Check your configuration."
    end if
    if (self%geo_class%provided .and. .not.self%exchange%geo_unit%required) then
      log_warn(*) "Input: geology class provided but not required. Check your configuration."
    end if
    if ((self%soil_class%provided .or. self%soil_horizon_class%provided) .and. .not.self%exchange%soil_id%required) then
      log_warn(*) "Input: soil class provided but not required. Check your configuration."
    end if
    if (self%lai_class%provided .and. .not.self%exchange%lai_class%required) then
      log_warn(*) "Input: LAI class provided but not required. Check your configuration."
    end if
    if (self%runoff%provided .and. .not.self%exchange%runoff_total%required) then
      log_warn(*) "Input: runoff provided but not required. Check your configuration."
    end if
    if (.not. self%exchange%runoff_total%required) self%runoff%provided = .false.
    call self%pre%reset_time(self%chunking, self%exchange%start_time)
    call self%pet%reset_time(self%chunking, self%exchange%start_time)
    call self%temp%reset_time(self%chunking, self%exchange%start_time)
    call self%tann%reset_time(self%chunking, self%exchange%start_time)
    call self%tmin%reset_time(self%chunking, self%exchange%start_time)
    call self%tmax%reset_time(self%chunking, self%exchange%start_time)
    call self%ssrd%reset_time(self%chunking, self%exchange%start_time)
    call self%strd%reset_time(self%chunking, self%exchange%start_time)
    call self%netrad%reset_time(self%chunking, self%exchange%start_time)
    call self%eabs%reset_time(self%chunking, self%exchange%start_time)
    call self%wind%reset_time(self%chunking, self%exchange%start_time)
    call self%runoff%reset_time(self%chunking, self%exchange%start_time)
  end subroutine input_initialize

  subroutine input_update(self)
    class(input_t), target, intent(inout) :: self
    log_trace(*) "Update Input"
    call self%pre%update(self%exchange%raw_pre%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%pet%update(self%exchange%raw_pet%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%temp%update(self%exchange%raw_temp%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%tann%update(self%exchange%raw_tann%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%tmin%update(self%exchange%raw_tmin%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%tmax%update(self%exchange%raw_tmax%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%ssrd%update(self%exchange%raw_ssrd%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%strd%update(self%exchange%raw_strd%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%netrad%update(self%exchange%raw_netrad%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%eabs%update(self%exchange%raw_eabs%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%wind%update(self%exchange%raw_wind%data, self%exchange%time, self%chunking, self%exchange%end_time)
    call self%runoff%update(self%exchange%runoff_total%data, self%exchange%time, self%chunking, self%exchange%end_time)
  end subroutine input_update

  subroutine input_finalize(self)
    class(input_t), target, intent(inout) :: self
    log_info(*) "Finalize Input"
    ! close datasets
    if (self%pre%provided .and. .not.self%pre%static) call self%pre%ds%close()
    if (self%pet%provided .and. .not.self%pet%static) call self%pet%ds%close()
    if (self%temp%provided .and. .not.self%temp%static) call self%temp%ds%close()
    if (self%tann%provided .and. .not.self%tann%static) call self%tann%ds%close()
    if (self%tmin%provided .and. .not.self%tmin%static) call self%tmin%ds%close()
    if (self%tmax%provided .and. .not.self%tmax%static) call self%tmax%ds%close()
    if (self%ssrd%provided .and. .not.self%ssrd%static) call self%ssrd%ds%close()
    if (self%strd%provided .and. .not.self%strd%static) call self%strd%ds%close()
    if (self%netrad%provided .and. .not.self%netrad%static) call self%netrad%ds%close()
    if (self%eabs%provided .and. .not.self%eabs%static) call self%eabs%ds%close()
    if (self%wind%provided .and. .not.self%wind%static) call self%wind%ds%close()
    if (self%runoff%provided) call self%runoff%ds%close()
    if (allocated(self%soil_class_one_layer)) deallocate(self%soil_class_one_layer)
    call self%exchange%lake_ids%clear(owned=.true.)
    call self%exchange%lake_max_levels%clear(owned=.true.)
    nullify(self%exchange%lake_points)
    self%lake_outlets = points_t()
    if (allocated(self%lake_ids)) deallocate(self%lake_ids)
    if (allocated(self%lake_max_levels)) deallocate(self%lake_max_levels)
    nullify(self%exchange%level0_lake)
    nullify(self%exchange%level0_land)
    if (self%owns_river_l0) then
      call self%exchange%fdir%clear(owned=.true.)
      if (associated(self%exchange%river_l0, self%river_l0)) nullify(self%exchange%river_l0)
      call self%river_l0%clean()
      self%owns_river_l0 = .false.
    end if
  end subroutine input_finalize

end module mo_input_container
