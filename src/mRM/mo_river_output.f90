!> \file    mo_river_output.f90
!> \copydoc mo_river_output

!> \brief   River output handler.
!> \details This module contains a methods to create NetCDF output of routed data.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    May 2025
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_output

  use mo_kind,   only: i1, i2, i4, i8, sp, dp
  use mo_constants, only : nodata_dp, nodata_sp, nodata_i4, nodata_i8
  use mo_river, only: river_t
  use mo_grid, only: cartesian
  use mo_grid_io, only: var, var_index, end_timestamp, time_values
  use mo_message, only: error_message, warn_message
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable, NcGroup, check
  use mo_datetime, only: datetime, timedelta, delta_from_string, decode_cf_time_units, one_day, one_hour
  use mo_utils, only: optval

  implicit none
  private
  public :: var, river_output_dataset

  !> \class river_output_variable
  !> \brief netcdf output variable container for a 2D variable
  type, extends(var) :: river_output_variable
    type(NcVariable) :: nc                     !< nc variable which contains the variable
    type(river_t), pointer :: river => null()  !< river the data is defined on (node based)
    real(sp), allocatable :: data_sp(:)        !< store the data between writes (real)
    real(dp), allocatable :: data_dp(:)        !< store the data between writes (real)
    integer(i4), allocatable :: data_i4(:)     !< store the data between writes (integer)
    integer(i8), allocatable :: data_i8(:)     !< store the data between writes (integer)
    integer(i4) :: counter = 0_i4              !< count the number of updateVariable calls
    logical :: static_written = .false.        !< static variable was written
  contains
    procedure, public :: init => out_var_init
    procedure, private :: out_var_update_sp, out_var_update_dp, out_var_update_i4, out_var_update_i8
    generic, public :: update => out_var_update_sp, out_var_update_dp, out_var_update_i4, out_var_update_i8
    procedure, public :: write => out_var_write
  end type river_output_variable

  !> \class river_output_dataset
  !> \brief netcdf output dataset handler for river based data
  !> \details Output dataset handler for static and temporal data on river nodes.
  type river_output_dataset
    type(river_t), pointer :: river => null()           !< river the data is defined on
    character(:), allocatable :: path                   !< path to the NetCDF file
    type(NcDataset) :: nc                               !< NcDataset to write
    type(river_output_variable), allocatable :: vars(:) !< store all created (dynamic) variables
    integer(i4) :: nvars                                !< number of variables
    logical :: static                                   !< only static variables (without time dimension)
    integer(i4) :: counter = 0_i4                       !< count written time steps
    type(datetime) :: previous_time                     !< previous time steps for bounds
    type(datetime) :: start_time                        !< start time for time units
    type(timedelta) :: delta                            !< time step in time units definition
    integer(i4) :: timestamp = end_timestamp            !< time stamp reference (0: begin, 1: center, 2: end of time interval)
    integer(i4) :: deflate_level = 6_i4                 !< deflate level for compression
  contains
    procedure, public :: init => output_init
    procedure, private :: output_update_sp, output_update_dp, output_update_i4, output_update_i8
    generic, public :: update => output_update_sp, output_update_dp, output_update_i4, output_update_i8
    procedure, public :: write => output_write
    procedure, public :: write_static => output_write_static
    procedure, public :: meta => output_meta
    procedure, public :: close => output_close
  end type river_output_dataset

contains

  !> \brief initialize river_output_variable
  subroutine out_var_init(self, meta, nc, river, dims, deflate_level)
    implicit none
    class(river_output_variable), intent(inout) :: self
    type(var), intent(in) :: meta !< variable definition
    type(NcDataset), intent(in) :: nc !< NcDataset to write
    type(river_t), pointer, intent(in) :: river !< river definition
    type(NcDimension), dimension(:), intent(in) :: dims !< dimensions in the file
    integer(i4), intent(in) :: deflate_level !< deflate level for compression

    self%name = meta%name
    self%static = meta%static
    self%avg = meta%avg
    self%dtype = "f64" ! default to double
    if (allocated(meta%dtype)) self%dtype = trim(meta%dtype)

    if (self%static) then
      self%nc = nc%setVariable(self%name, self%dtype, dims(:1), deflate_level=deflate_level, shuffle=.true.)
    else
      self%nc = nc%setVariable(self%name, self%dtype, dims(:2), deflate_level=deflate_level, shuffle=.true.)
    end if

    if (allocated(meta%long_name)) self%long_name = meta%long_name
    if (allocated(meta%standard_name)) self%standard_name = meta%standard_name
    if (allocated(meta%units)) self%units = meta%units

    if (allocated(self%long_name)) call self%nc%setAttribute("long_name", self%long_name)
    if (allocated(self%standard_name)) call self%nc%setAttribute("standard_name", self%standard_name)
    if (allocated(self%units)) call self%nc%setAttribute("units", self%units)

    ! set coordinates
    call self%nc%setAttribute("coordinates", "x y")
    call self%nc%setAttribute("mesh", "river")
    call self%nc%setAttribute("location", "node")

    self%river => river
    ! input data is still either real(dp) or integer(i4)
    select case(self%dtype)
      case("f32")
        call self%nc%setFillValue(nodata_sp)
        call self%nc%setAttribute("missing_value", nodata_sp)
        self%kind = "dp" ! double precision as default for real data
      case("f64")
        call self%nc%setFillValue(nodata_dp)
        call self%nc%setAttribute("missing_value", nodata_dp)
        self%kind = "dp" ! double precision as default for real data
      case("i32")
        call self%nc%setFillValue(nodata_i4)
        call self%nc%setAttribute("missing_value", nodata_i4)
        self%kind = "i4" ! 4 byte integers as default for int data
      case("i64")
        call self%nc%setFillValue(nodata_i8)
        call self%nc%setAttribute("missing_value", nodata_i8)
        self%kind = "i4" ! 4 byte integers as default for int data
      case default
        call error_message("river_output_variable: unsupported dtype: ", self%name, ": ", self%dtype)
    end select
    if (allocated(meta%kind)) then
      self%kind = meta%kind
      if ((self%dtype(1:1) == "f" .and. meta%kind(2:2) /= "p") .or. (self%dtype(1:1) == "i" .and. meta%kind(1:1) /= "i")) &
        call warn_message("river_output_variable: variable dtype and array kind will result in conversion: ", &
                          self%name, ", dtype: ", self%dtype, ", kind:", self%kind)
    end if
    select case(self%kind)
      case("sp")
        allocate(self%data_sp(self%river%n_nodes), source=0.0_sp)
      case("dp")
        allocate(self%data_dp(self%river%n_nodes), source=0.0_dp)
      case("i4")
        allocate(self%data_i4(self%river%n_nodes), source=0_i4)
      case("i8")
        allocate(self%data_i8(self%river%n_nodes), source=0_i8)
      case default
        call error_message("river_output_variable: unsupported kind: ", self%name, ": ", self%kind)
    end select
  end subroutine out_var_init

  !> \brief Update river_output_variable
  subroutine out_var_update_sp(self, data)
    implicit none
    class(river_output_variable), intent(inout) :: self
    real(sp), intent(in), dimension(:) :: data !< data for current time step
    if (.not.allocated(self%data_sp)) call error_message("river_output_variable: wrong kind: ", self%name, ", ", self%kind, "=/=sp")
    self%data_sp = self%data_sp + data
    self%counter = self%counter + 1_i4
  end subroutine out_var_update_sp

  !> \brief Update river_output_variable
  subroutine out_var_update_dp(self, data)
    implicit none
    class(river_output_variable), intent(inout) :: self
    real(dp), intent(in), dimension(:) :: data !< data for current time step
    if (.not.allocated(self%data_dp)) call error_message("river_output_variable: wrong kind: ", self%name, ", ", self%kind, "=/=dp")
    self%data_dp = self%data_dp + data
    self%counter = self%counter + 1_i4
  end subroutine out_var_update_dp

  !> \brief Update river_output_variable
  subroutine out_var_update_i4(self, data)
    implicit none
    class(river_output_variable), intent(inout) :: self
    integer(i4), intent(in), dimension(:) :: data !< data for current time step
    if (.not.allocated(self%data_i4)) call error_message("river_output_variable: wrong kind: ", self%name, ", ", self%kind, "=/=i4")
    self%data_i4 = self%data_i4 + data
    self%counter = self%counter + 1_i4
  end subroutine out_var_update_i4

  !> \brief Update river_output_variable
  subroutine out_var_update_i8(self, data)
    implicit none
    class(river_output_variable), intent(inout) :: self
    integer(i8), intent(in), dimension(:) :: data !< data for current time step
    if (.not.allocated(self%data_i8)) call error_message("river_output_variable: wrong kind: ", self%name, ", ", self%kind, "=/=i8")
    self%data_i8 = self%data_i8 + data
    self%counter = self%counter + 1_i4
  end subroutine out_var_update_i8

  !> \brief Write timestep to file
  !> \details Write the content of the derived types's component 'data' to file, average if necessary
  !> \changelog
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !!
  !> \authors David Schafer
  !> \date June 2015
  subroutine out_var_write(self, t_index)
    implicit none
    class(river_output_variable), intent(inout) :: self
    !> index along the time dimension of the netcdf variable
    integer(i4), intent(in), optional :: t_index
    if (self%counter == 0_i4) call error_message("river_output_variable: no data was added before writing: ", self%name)
    if (self%static .and. self%static_written) return
    if (self%avg) then
      select case(self%kind)
        case("sp")
          self%data_sp = self%data_sp / real(self%counter, sp)
        case("dp")
          self%data_dp = self%data_dp / real(self%counter, dp)
        case("i4")
          self%data_i4 = self%data_i4 / self%counter ! this is rounding
        case("i8")
          self%data_i8 = self%data_i8 / int(self%counter, i8)  ! this is rounding
      end select
    end if
    if (self%static) then
      select case(self%kind)
        case("sp")
          call self%nc%setData(self%data_sp)
        case("dp")
          call self%nc%setData(self%data_dp)
        case("i4")
          call self%nc%setData(self%data_i4)
        case("i8")
          call self%nc%setData(self%data_i8)
      end select
      self%static_written = .true.
    else
      if (.not.present(t_index)) call error_message("river_output_variable: no time index was given for temporal variable: ", self%name)
      select case(self%kind)
        case("sp")
          call self%nc%setData(self%data_sp, [1_i4, t_index])
        case("dp")
          call self%nc%setData(self%data_dp, [1_i4, t_index])
        case("i4")
          call self%nc%setData(self%data_i4, [1_i4, t_index])
        case("i8")
          call self%nc%setData(self%data_i8, [1_i4, t_index])
      end select
    end if
    ! reset
    select case(self%kind)
      case("sp")
        self%data_sp = 0.0_sp
      case("dp")
        self%data_dp = 0.0_dp
      case("i4")
        self%data_i4 = 0_i4
      case("i8")
        self%data_i8 = 0_i8
    end select
    self%counter = 0_i4
  end subroutine out_var_write

  !> \brief Initialize river_output_dataset
  !> \details Create and initialize the output file handler.
  !> \authors Matthias Zink
  !> \authors Robert Schweppe
  !> \authors Sebastian Müller
  !> \date Apr 2013
  subroutine output_init(self, path, river, vars, start_time, delta, timestamp, deflate_level, network)
    implicit none
    class(river_output_dataset), intent(inout) :: self
    character(*), intent(in) :: path !< path to the file
    type(river_t), pointer, intent(in) :: river !< river definition for this output file
    type(var), dimension(:), intent(in) :: vars !< variables of the output file
    type(datetime), intent(in), optional :: start_time !< reference time
    character(*), intent(in), optional :: delta !< time units delta ("minutes", "hours" (default), "days")
    integer(i4), intent(in), optional :: timestamp !< time stamp location in time span (0: begin, 1: center, 2: end (default))
    integer(i4), intent(in), optional :: deflate_level !< deflate level for compression
    logical, intent(in), optional :: network !< store river-network in the file (default: .true.)

    character(:), allocatable :: units, units_dt
    type(NcDimension) :: t_dim, two_dim, node_dim, link_dim, dims(2)
    type(NcVariable) :: node_x_var, node_y_var, t_var, mesh_var, link_var
    integer(i4) :: i, nlinks, topo_dim
    integer(i8) :: n, k
    integer(i8), allocatable :: links(:, :)
    logical :: net

    net = optval(network, .true.)
    topo_dim = merge(1_i4, 0_i4, net)

    self%path = trim(path)
    self%nc = NcDataset(self%path, "w")
    self%river => river
    self%counter = 0_i4

    self%deflate_level = optval(deflate_level, 6_i4)
    self%timestamp = optval(timestamp, end_timestamp)

    self%nvars = size(vars)
    self%static = .true.
    do i = 1_i4, self%nvars
      self%static = self%static .and. vars(i)%static
    end do

    ! node_dim = self%nc%setDimension("node") ! use unlimited dimension to support i8 index
    ! link_dim = self%nc%setDimension("link") ! use unlimited dimension to support i8 index
    two_dim = self%nc%setDimension("Two", 2_i4)
    node_dim = self%nc%setDimension("node", int(self%river%n_nodes, i4))  ! only works if network is not to huge for i4

    ! 1D network topology following UGRID conventions
    mesh_var = self%nc%setVariable("river", "i32", dims(:0)) ! mesh variable as scalar integer
    call mesh_var%setAttribute("cf_role", "mesh_topology")
    call mesh_var%setAttribute("long_name", "river network definition")
    call mesh_var%setAttribute("topology_dimension", topo_dim)  ! 0 - only nodes, 1 - with links
    call mesh_var%setAttribute("node_coordinates", "river_node_x river_node_y")
    if (net) call mesh_var%setAttribute("edge_node_connectivity", "links")

    ! coordinates
    node_x_var = self%nc%setVariable("river_node_x", "f64", [node_dim])
    node_y_var = self%nc%setVariable("river_node_y", "f64", [node_dim])

    if (self%river%grid%coordsys==cartesian) then
      call node_x_var%setAttribute("long_name", "x coordinate of projection")
      call node_x_var%setAttribute("standard_name", "projection_x_coordinate")
      call node_x_var%setAttribute("units", "m") ! TODO: this should be configurable
      call node_y_var%setAttribute("long_name", "y coordinate of projection")
      call node_y_var%setAttribute("standard_name", "projection_y_coordinate")
      call node_y_var%setAttribute("units", "m") ! TODO: this should be configurable
    else
      call node_x_var%setAttribute("long_name", "longitude")
      call node_x_var%setAttribute("standard_name", "longitude")
      call node_x_var%setAttribute("units", "degrees_east")
      call node_y_var%setAttribute("long_name", "latitude")
      call node_y_var%setAttribute("standard_name", "latitude")
      call node_y_var%setAttribute("units", "degrees_north")
    end if
    ! this should set the node-dim size
    call node_x_var%setData(self%river%points%x)
    call node_y_var%setData(self%river%points%y)

    ! links
    if (net) then
      nlinks = int(self%river%n_nodes, i4) - size(self%river%sinks, kind=i4)
      link_dim = self%nc%setDimension("link", nlinks)
      link_var = self%nc%setVariable("links", "i64", [two_dim, link_dim])
      call link_var%setAttribute("cf_role", "edge_node_connectivity")
      call link_var%setAttribute("long_name", "river links definition")
      call link_var%setAttribute("start_index", 1)  ! fortran indices starting with 1
      allocate(links(2_i4, nlinks))
      k = 0_i8
      do n = 1_i8, self%river%n_nodes
        if (self%river%is_sink(n)) cycle
        k = k + 1_i8
        links(1_i8, k) = n
        links(2_i8, k) = self%river%down(n)
      end do
      call link_var%setData(links)
      deallocate(links)
    end if

    if (.not.self%static) then
      if (.not.present(start_time)) call error_message("output: if dataset is not static, a start_time is needed")
      units_dt = trim(optval(delta, "hours"))
      self%previous_time = start_time
      self%start_time = start_time
      self%delta = delta_from_string(units_dt)
      units = units_dt // " since " // start_time%str()
      t_dim = self%nc%setDimension("time", 0)
      t_var = self%nc%setVariable("time", "i32", [t_dim])
      call t_var%setAttribute("long_name", "time")
      call t_var%setAttribute("standard_name", "time")
      call t_var%setAttribute("axis", "T")
      call t_var%setAttribute("units", units)
      call t_var%setAttribute("bounds", "time_bnds")
      t_var = self%nc%setVariable("time_bnds", "i32", [two_dim, t_dim])
    end if

    dims(:) = [node_dim, t_dim]

    allocate(self%vars(self%nvars))
    do i = 1_i4, self%nvars
      call self%vars(i)%init(vars(i), self%nc, self%river, dims, self%deflate_level)
    end do
  end subroutine output_init

  !> \brief Update a variable
  !> \details Add the array given as actual argument to the derived type's component 'data'
  subroutine output_update_sp(self, name, data)
    implicit none
    class(river_output_dataset), intent(inout) :: self
    character(*), intent(in) :: name !< name of the variable
    real(sp), intent(in), dimension(:) :: data !< data for current time step
    call self%vars(var_index(self%vars, name, "output%update"))%update(data)
  end subroutine output_update_sp

  !> \brief Update a variable
  !> \details Add the array given as actual argument to the derived type's component 'data'
  subroutine output_update_dp(self, name, data)
    implicit none
    class(river_output_dataset), intent(inout) :: self
    character(*), intent(in) :: name !< name of the variable
    real(dp), intent(in), dimension(:) :: data !< data for current time step
    call self%vars(var_index(self%vars, name, "output%update"))%update(data)
  end subroutine output_update_dp

  !> \brief Update a variable
  !> \details Add the array given as actual argument to the derived type's component 'data'
  subroutine output_update_i4(self, name, data)
    implicit none
    class(river_output_dataset), intent(inout) :: self
    character(*), intent(in) :: name !< name of the variable
    integer(i4), intent(in), dimension(:) :: data !< data for current time step
    call self%vars(var_index(self%vars, name, "output%update"))%update(data)
  end subroutine output_update_i4

  !> \brief Update a variable
  !> \details Add the array given as actual argument to the derived type's component 'data'
  subroutine output_update_i8(self, name, data)
    implicit none
    class(river_output_dataset), intent(inout) :: self
    character(*), intent(in) :: name !< name of the variable
    integer(i8), intent(in), dimension(:) :: data !< data for current time step
    call self%vars(var_index(self%vars, name, "output%update"))%update(data)
  end subroutine output_update_i8

  !> \brief Write all accumulated data.
  !> \details Write all accumulated and potentially averaged data to disk.
  !> \changelog
  !! - Robert Schweppe Jun 2018
  !!   - refactoring and reformatting
  !!
  !> \authors David Schaefer
  !> \date June 2015
  subroutine output_write(self, current_time)
    implicit none
    class(river_output_dataset), intent(inout) :: self
    type(datetime), intent(in), optional :: current_time !< end time of the current time span
    integer(i4) :: i, t_start, t_end, t_stamp
    type(NcVariable) :: t_var
    self%counter = self%counter + 1_i4
    ! add to time variable
    if (.not.self%static) then
      if (.not.present(current_time)) call error_message("output: no time was given: ", self%path)
      call time_values(self%start_time, self%previous_time, current_time, self%delta, self%timestamp, &
                       t_start, t_end, t_stamp) ! intent(out)
      t_var = self%nc%getVariable("time")
      call t_var%setData(t_stamp, [self%counter])
      t_var = self%nc%getVariable("time_bnds")
      call t_var%setData(t_start, [1, self%counter])
      call t_var%setData(t_end, [2, self%counter])
      self%previous_time = current_time
    end if
    ! write all variables
    do i = 1_i4, self%nvars
      call self%vars(i)%write(self%counter)
    end do
  end subroutine output_write

  !> \brief Write all accumulated static data.
  subroutine output_write_static(self)
    implicit none
    class(river_output_dataset), intent(inout) :: self
    integer(i4) :: i
    do i = 1_i4, self%nvars
      if (self%vars(i)%static) call self%vars(i)%write()
    end do
  end subroutine output_write_static

  !> \brief Get variable meta data.
  !> \return \ref var meta data definition
  type(var) function output_meta(self, name)
    implicit none
    class(river_output_dataset) :: self
    character(*), intent(in) :: name !< name of the variable
    output_meta = self%vars(var_index(self%vars, name, "output%meta"))%meta()
  end function output_meta

  !> \brief Close the file
  subroutine output_close(self)
    implicit none
    class(river_output_dataset) :: self
    integer(i4) :: i
    do i = 1_i4, self%nvars
      ! check if variables have buffered data that was not written
      if (self%vars(i)%counter > 0_i4) call warn_message("output%close: unwritten buffered data for variable: ", self%vars(i)%name)
    end do
    call self%nc%close()
    deallocate(self%vars)
  end subroutine output_close

end module mo_river_output
