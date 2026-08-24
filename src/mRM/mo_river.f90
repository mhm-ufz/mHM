!> \file    mo_river.f90
!> \copydoc mo_river

!> \brief   River representation.
!> \details This module contains a river network representation based on DAG and the grid class.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    May 2025
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river

  use mo_kind,   only: i1, i2, i4, i8, sp, dp
  use mo_dag, only: branching, order_t, traversal_visit
  use mo_grid, only: grid_t, bottom_up, cartesian, dist_latlon
  use mo_grid_io, only: var, output_dataset
  use mo_message, only: error_message, message
  use mo_utils, only: optval
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable
  use mo_constants, only: nodata_i1, nodata_i2, nodata_i4, nodata_i8, nodata_dp

  implicit none
  private

  !> \name D8 direction values
  !> \brief Constants describing the D8 routing direction of a cell.
  !!@{
  integer(i2), public, parameter :: sink = 0_i2 !< sink
  integer(i2), public, parameter :: d8_E = 1_i2 !< east
  integer(i2), public, parameter :: d8_SE = 2_i2 !< south-east
  integer(i2), public, parameter :: d8_S = 4_i2 !< south
  integer(i2), public, parameter :: d8_SW = 8_i2 !< south-west
  integer(i2), public, parameter :: d8_W = 16_i2 !< west
  integer(i2), public, parameter :: d8_NW = 32_i2 !< north-west
  integer(i2), public, parameter :: d8_N = 64_i2 !< north
  integer(i2), public, parameter :: d8_NE = 128_i2 !< north-east
  !> all directions as array
  integer(i2), public, dimension(8), parameter :: d8_all = [d8_E, d8_SE, d8_S, d8_SW, d8_W, d8_NW, d8_N, d8_NE]
  !> matching back pointing directions as array
  integer(i2), public, dimension(8), parameter :: d8_back = [d8_W, d8_NW, d8_N, d8_NE, d8_E, d8_SE, d8_S, d8_SW]
  !> leaving directions at east border as array
  integer(i2), public, dimension(3), parameter :: d8_leave_E = [d8_NE, d8_E, d8_SE]
  !> leaving directions at south border as array
  integer(i2), public, dimension(3), parameter :: d8_leave_S = [d8_SE, d8_S, d8_SW]
  !> leaving directions at west border as array
  integer(i2), public, dimension(3), parameter :: d8_leave_W = [d8_SW, d8_W, d8_NW]
  !> leaving directions at north border as array
  integer(i2), public, dimension(3), parameter :: d8_leave_N = [d8_NW, d8_N, d8_NE]
  !> leaving directions at south-east corner as array
  integer(i2), public, dimension(5), parameter :: d8_leave_SE = [d8_NE, d8_E, d8_SE, d8_S, d8_SW]
  !> leaving directions at south-west corner as array
  integer(i2), public, dimension(5), parameter :: d8_leave_SW = [d8_SE, d8_S, d8_SW, d8_W, d8_NW]
  !> leaving directions at north-east corner as array
  integer(i2), public, dimension(5), parameter :: d8_leave_NE = [d8_NW, d8_N, d8_NE, d8_E, d8_SE]
  !> leaving directions at north-west corner as array
  integer(i2), public, dimension(5), parameter :: d8_leave_NW = [d8_SW, d8_W, d8_NW, d8_N, d8_NE]
  !!@}

  !> \name ldd direction values
  !> \brief Constants describing the Local drain direction (lld) routing direction of a cell.
  !!@{
  integer(i1), public, parameter :: o = 0_i1 !< no-data shortcut
  integer(i1), public, parameter :: ldd_sink = 5_i1 !< sink
  integer(i1), public, parameter :: ldd_E = 6_i1 !< east
  integer(i1), public, parameter :: ldd_SE = 3_i1 !< south-east
  integer(i1), public, parameter :: ldd_S = 2_i1 !< south
  integer(i1), public, parameter :: ldd_SW = 1_i1 !< south-west
  integer(i1), public, parameter :: ldd_W = 4_i1 !< west
  integer(i1), public, parameter :: ldd_NW = 7_i1 !< north-west
  integer(i1), public, parameter :: ldd_N = 8_i1 !< north
  integer(i1), public, parameter :: ldd_NE = 9_i1 !< north-east
  !> all directions as array
  integer(i1), public, dimension(8), parameter :: ldd_all = [ldd_E, ldd_SE, ldd_S, ldd_SW, ldd_W, ldd_NW, ldd_N, ldd_NE]
  !> matching back pointing directions as array
  integer(i1), public, dimension(8), parameter :: ldd_back = [ldd_W, ldd_NW, ldd_N, ldd_NE, ldd_E, ldd_SE, ldd_S, ldd_SW]
  !> x direction change for given ldd flow direction
  integer(i4), public, dimension(9), parameter :: ldd_dx = [-1_i4, 0_i4, 1_i4, -1_i4, 0_i4, 1_i4, -1_i4, 0_i4, 1_i4]
  !> y direction change for given ldd flow direction (bottom-up)
  integer(i4), public, dimension(9), parameter :: ldd_dy = [-1_i4, -1_i4, -1_i4, 0_i4, 0_i4, 0_i4, 1_i4, 1_i4, 1_i4]
  !> ldd to d8 conversion
  integer(i2), public, dimension(9), parameter :: map_ldd_to_d8 = [d8_SW, d8_S, d8_SE, d8_W, sink, d8_E, d8_NW, d8_N, d8_NE]
  !> d8 to ldd conversion (array mostly empty)
  integer(i1), public, dimension(128), parameter :: map_d8_to_ldd = [ &
    ldd_E, &
    ldd_SE, o, &
    ldd_S,  o,o,o, &
    ldd_SW, o,o,o,o,o,o,o, &
    ldd_W,  o,o,o,o,o,o,o,o,o,o,o,o,o,o,o, &
    ldd_NW, o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o, &
    ldd_N,  o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o, &
            o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o,o, &
    ldd_NE]
  !!@}

  !> \class river_t
  !> \brief River network representation
  type, extends(branching), public :: river_t
    type(grid_t), pointer :: grid => null() !< grid the river network is defined on
    logical, allocatable :: is_sink(:) !< flag to indicate sinks size(n_nodes)
    integer(i2), allocatable :: fdir(:) !< D8 flow direction (only for a D8-river) size(ncells)
    integer(i4), allocatable :: facc(:) !< flow accumulation size(n_nodes)
    real(dp), allocatable :: upstream_area(:) !< upstream area of node size(n_nodes)
    real(dp), allocatable :: link_length(:) !< length of link starting at node (0 if node is sink) size(n_nodes)
    real(dp), allocatable :: link_slope(:) !< slope of link starting at node (in %) size(n_nodes)
    real(dp), allocatable :: celerity(:) !< celerity of link starting at node set by upscaler size(n_nodes)
    type(order_t) :: order !< level based order of the network
    ! scc related attributes
    logical :: scc = .false. !< indicate that this river is a SCC-river (not D8)
    integer(i8), allocatable :: node_cell(:) !< map node to grid cell id size(n_nodes)
    integer(i8), allocatable :: cell_node_select(:) !< select contained node to represent cell (e.g. highest facc) size(ncells)
    real(dp), allocatable :: area_fraction(:) !< area fraction for each node in cell size(n_nodes)
    real(dp), allocatable :: node_x(:) !< x coordinate for each node size(n_nodes)
    real(dp), allocatable :: node_y(:) !< x coordinate for each node size(n_nodes)
  contains
    procedure, public :: from_fdir => river_from_fdir
    procedure, public :: calc_order => river_order
    procedure, public :: calc_fdir => river_fdir
    procedure, public :: calc_facc => river_facc
    procedure, public :: label_subcatchments => river_label_subcatchments
    procedure, public :: calc_upstream_area => river_upstream_area
    procedure, public :: calc_length => river_length
    procedure, public :: calc_slope => river_slope
    procedure, public :: calc_celerity => river_celerity
    procedure, public :: select_cell_values => river_select_cell_values
    procedure, public :: export => river_export
    procedure, public :: clean => river_destroy
    procedure, public :: to_restart_dataset => river_to_restart_dataset, to_restart_file => river_to_restart_file
    generic, public :: to_restart => to_restart_dataset, to_restart_file
    procedure, public :: from_restart_dataset => river_from_restart_dataset, from_restart_file => river_from_restart_file
    generic, public :: from_restart => from_restart_dataset, from_restart_file
  end type river_t

contains

  !> \brief Convert D8 flow direction to LDD flow direction
  subroutine d8_to_ldd(fdir_d8, fdir_ldd)
    integer(i2), dimension(:), intent(in) :: fdir_d8 !< D8 flow direction
    integer(i1), dimension(:), intent(out) :: fdir_ldd !< LDD flow direction
    integer(i8) :: i
    if (size(fdir_ldd, kind=i8) /= size(fdir_d8, kind=i8)) then
      call error_message("d8_to_ldd: input and output arrays have different size")
    end if
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i8, size(fdir_d8, kind=i8)
      select case (fdir_d8(i))
        case (sink)
          fdir_ldd(i) = ldd_sink
        case (d8_E)
          fdir_ldd(i) = ldd_E
        case (d8_SE)
          fdir_ldd(i) = ldd_SE
        case (d8_S)
          fdir_ldd(i) = ldd_S
        case (d8_SW)
          fdir_ldd(i) = ldd_SW
        case (d8_W)
          fdir_ldd(i) = ldd_W
        case (d8_NW)
          fdir_ldd(i) = ldd_NW
        case (d8_N)
          fdir_ldd(i) = ldd_N
        case (d8_NE)
          fdir_ldd(i) = ldd_NE
        case default
          fdir_ldd(i) = nodata_i1
      end select
    end do
    !$omp end parallel do
  end subroutine d8_to_ldd

  !> \brief Convert LDD flow direction to D8 flow direction
  subroutine ldd_to_d8(fdir_ldd, fdir_d8)
    integer(i1), dimension(:), intent(in) :: fdir_ldd !< LDD flow direction
    integer(i2), dimension(:), intent(out) :: fdir_d8 !< D8 flow direction
    integer(i8) :: i
    if (size(fdir_ldd, kind=i8) /= size(fdir_d8, kind=i8)) then
      call error_message("ldd_to_d8: input and output arrays have different size")
    end if
    !$omp parallel do default(shared) schedule(static) private(i)
    do i = 1_i8, size(fdir_d8, kind=i8)
      select case (fdir_ldd(i))
        case (ldd_sink)
          fdir_d8(i) = sink
        case (ldd_E)
          fdir_d8(i) = d8_E
        case (ldd_SE)
          fdir_d8(i) = d8_SE
        case (ldd_S)
          fdir_d8(i) = d8_S
        case (ldd_SW)
          fdir_d8(i) = d8_SW
        case (ldd_W)
          fdir_d8(i) = d8_W
        case (ldd_NW)
          fdir_d8(i) = d8_NW
        case (ldd_N)
          fdir_d8(i) = d8_N
        case (ldd_NE)
          fdir_d8(i) = d8_NE
        case default
          fdir_d8(i) = nodata_i2
      end select
    end do
    !$omp end parallel do
  end subroutine ldd_to_d8

  !> \brief Get all links to and from given cell depending on flow direction
  pure subroutine get_links(mask, cells, fdir, i, j, dy, periodic, n_up, up, down)
    logical, dimension(:,:), intent(in) :: mask !< grid mask
    integer(i8), dimension(:,:), intent(in) :: cells !< cells id matrix
    integer(i2), dimension(:), intent(in) :: fdir !< flow direction array
    integer(i4), intent(in) :: i !< i index of current cell (on x-axis)
    integer(i4), intent(in) :: j !< j index of current cell (on y-axis)
    integer(i4), intent(in) :: dy !< direction of north in the grid matrix (1/-1)
    logical, intent(in) :: periodic !< flag to indicate periodic latlon grid (360 deg lon axis)
    integer(i4), intent(out) :: n_up !< number of upstream neighbor cells
    integer(i8), dimension(8), intent(out) :: up !< all upstream neighbor cells
    integer(i8), intent(out) :: down !< the downstream cell
    integer(i4) :: imax, jmax, ni, nj, k
    imax = size(cells, dim=1)
    jmax = size(cells, dim=2)
    ! upstream
    n_up = 0_i4
    do k = 1_i4, 8_i4
      call next(d8_all(k), dy, i, j, ni, nj) ! check all directions
      if (periodic) then ! sinks still indicated by nj=0
        if (ni < 1_i4) ni = imax
        if (imax < ni) ni = 1_i4
      end if
      if (ni < 1_i4 .or. imax < ni .or. nj < 1_i4 .or. jmax < nj ) cycle ! outside matrix / sink
      if (.not.mask(ni,nj)) cycle ! outside mask
      if (fdir(cells(ni,nj)) /= d8_back(k)) cycle ! check if this next cell is pointing back
      n_up = n_up + 1_i4 ! found an inflow
      up(n_up) = cells(ni,nj) ! add the upstream neighbor to the list
    end do
    ! downstream
    call get_down(mask, cells, fdir(cells(i,j)), i, j, dy, periodic, down)
  end subroutine get_links

  pure subroutine get_down(mask, cells, fdir, i, j, dy, periodic, down)
    logical, dimension(:,:), intent(in) :: mask !< grid mask
    integer(i8), dimension(:,:), intent(in) :: cells !< cells id matrix
    integer(i2), intent(in) :: fdir !< flow direction of current cell
    integer(i4), intent(in) :: i !< i index of current cell (on x-axis)
    integer(i4), intent(in) :: j !< j index of current cell (on y-axis)
    integer(i4), intent(in) :: dy !< direction of north in the grid matrix (1/-1)
    logical, intent(in) :: periodic !< flag to indicate periodic latlon grid (360 deg lon axis)
    integer(i8), intent(out) :: down !< the downstream cell
    integer(i4) :: imax, jmax, ni, nj
    imax = size(cells, dim=1)
    jmax = size(cells, dim=2)
    down = 0_i8
    call next(fdir, dy, i, j, ni, nj) ! get downstream cell
    if (periodic) then ! sinks still indicated by nj=0
      if (ni < 1_i4) ni = imax
      if (imax < ni) ni = 1_i4
    end if
    if (ni < 1_i4 .or. imax < ni .or. nj < 1_i4 .or. jmax < nj ) return ! outside matrix / sink
    if (.not.mask(ni,nj)) return ! outside mask
    down = cells(ni,nj)
  end subroutine get_down

  !> \brief Get next matrix indices from flow direction: (i,j) -> (ni,nj)
  pure subroutine next(fdir, dy, i, j, ni, nj)
    integer(i2), intent(in) :: fdir !< flow direction
    integer(i4), intent(in) :: dy !< direction of north in the grid matrix (1/-1)
    integer(i4), intent(in) :: i !< i index of current cell (on x-axis)
    integer(i4), intent(in) :: j !< j index of current cell (on y-axis)
    integer(i4), intent(out) :: ni !< i index of next cell (on x-axis)
    integer(i4), intent(out) :: nj !< j index of next cell (on y-axis)
    select case(fdir)
      case(d8_E) ! E
        ni = i+1_i4
        nj = j
      case(d8_SE) ! SE
        ni = i+1_i4
        nj = j-dy
      case(d8_S) ! S
        ni = i
        nj = j-dy
      case(d8_SW) ! SW
        ni = i-1_i4
        nj = j-dy
      case(d8_W) ! W
        ni = i-1_i4
        nj = j
      case(d8_NW) ! NW
        ni = i-1_i4
        nj = j+dy
      case(d8_N) ! N
        ni = i
        nj = j+dy
      case(d8_NE) ! NE
        ni = i+1_i4
        nj = j+dy
      case default ! sink/outlet
        ni = 0_i4
        nj = 0_i4
    end select
  end subroutine next

  !> \brief Get next matrix indices from flow direction: (i,j) -> (ni,nj)
  pure integer(i2) function get_fdir(from, to, dy)
    integer(i4), intent(in) :: from(2) !< from-cell indices
    integer(i4), intent(in) :: to(2) !< to-cell indices
    integer(i4), intent(in) :: dy !< direction of north in the grid matrix (1/-1)
    if (to(1)==(from(1)+1_i4) .and. to(2)==from(2)) then
      get_fdir = d8_E
    else if (to(1)==(from(1)+1_i4) .and. to(2)==(from(2)-dy)) then
      get_fdir = d8_SE
    else if (to(1)==from(1) .and. to(2)==(from(2)-dy)) then
      get_fdir = d8_S
    else if (to(1)==(from(1)-1_i4) .and. to(2)==(from(2)-dy)) then
      get_fdir = d8_SW
    else if (to(1)==(from(1)-1_i4) .and. to(2)==from(2)) then
      get_fdir = d8_W
    else if (to(1)==(from(1)-1_i4) .and. to(2)==(from(2)+dy)) then
      get_fdir = d8_NW
    else if (to(1)==from(1) .and. to(2)==(from(2)+dy)) then
      get_fdir = d8_N
    else if (to(1)==(from(1)+1_i4) .and. to(2)==(from(2)+dy)) then
      get_fdir = d8_NE
    else
      get_fdir = sink
    end if
  end function get_fdir

  !> \brief Initialize river network from flow direction.
  subroutine river_from_fdir(this, fdir, grid, calculate_length, calculate_node_xy)
    class(river_t), intent(inout) :: this
    integer(i2), dimension(:), intent(in) :: fdir !< D8 flow direction
    type(grid_t), pointer, intent(in), optional :: grid !< grid the river network is defined on
    logical, intent(in), optional :: calculate_length !< whether to calculate the link length from fdir (default: .true.)
    logical, intent(in), optional :: calculate_node_xy !< whether to calculate node locations (default: .true.)
    integer(i8), allocatable :: cells(:,:), down(:)
    integer(i8) :: downstream,  i
    integer(i4) :: dy ! direction of north in the grid matrix (1/-1)
    logical :: periodic, calc_length, calc_node_xy
    real(dp), allocatable :: xax(:), yax(:)

    calc_length = optval(calculate_length, .true.)
    calc_node_xy = optval(calculate_node_xy, .true.)
    if (present(grid)) this%grid => grid
    if (.not.associated(this%grid)) call error_message("river%from_fdir: grid not associated")
    if (size(fdir, kind=i8) /= this%grid%ncells) call error_message("river%from_fdir: size of fdir does not match grid ncells")

    this%scc = .false. ! if derived from fdir, this is a D8-river
    allocate(this%fdir(this%grid%ncells))
    ! allocate(this%node_cell(this%grid%ncells)) ! nodes correspond to cells
    ! allocate(this%cell_node_select(this%grid%ncells)) ! nodes correspond to cells
    ! allocate(this%area_fraction(this%grid%ncells)) ! D8 has no area fraction
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i8, this%grid%ncells
      this%fdir(i) = fdir(i)
      ! this%node_cell(i) = i
      ! this%cell_node_select(i) = i
      ! this%area_fraction(i) = 1.0_dp
    end do
    !$omp end parallel do

    periodic = this%grid%is_periodic()
    dy = -1_i4 ! top-down grid starts north
    if (this%grid%y_direction==bottom_up) dy = 1_i4

    allocate(cells(this%grid%nx, this%grid%ny))
    call this%grid%gen_id_matrix(cells)
    allocate(down(this%grid%ncells))

    !$omp parallel do default(shared) private(downstream) schedule(static)
    do i = 1_i8, this%grid%ncells
      call get_down(this%grid%mask, cells, this%fdir(i), &
                    this%grid%cell_ij(i,1), this%grid%cell_ij(i,2), dy, periodic, &
                    downstream)
      down(i) = downstream
      if (downstream == 0_i8 .and. this%fdir(i) /= sink) then
        this%fdir(i) = sink ! ensure sink direction
      end if
    end do
    !$omp end parallel do
    deallocate(cells)
    call this%init(down)
    deallocate(down)

    ! determine sinks
    allocate(this%is_sink(this%grid%ncells))
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i8, this%grid%ncells
      this%is_sink(i) = (this%down(i) == 0_i8)
    end do
    !$omp end parallel do

    ! calculate d8 reach length
    if (calc_length) call this%calc_length()

    ! determine node locations
    if (calc_node_xy) then
      allocate(this%node_x(this%n_nodes))
      allocate(this%node_y(this%n_nodes))
      xax = this%grid%x_axis()
      yax = this%grid%y_axis()
      !$omp parallel do default(shared) schedule(static)
      do i = 1_i8, this%grid%ncells
        this%node_x(i) = xax(this%grid%cell_ij(i, 1))
        this%node_y(i) = yax(this%grid%cell_ij(i, 2))
      end do
      !$omp end parallel do
    end if
  end subroutine river_from_fdir

  !> \brief Get river node order by levels
  subroutine river_order(this, root)
    class(river_t), intent(inout) :: this
    logical, optional, intent(in) :: root !< level as distance from root (.true.) or max. distance from headwater (.false., default)
    integer(i8) :: istat
    if (allocated(this%order%id)) deallocate(this%order%id)
    if (allocated(this%order%level_start)) deallocate(this%order%level_start)
    if (allocated(this%order%level_end)) deallocate(this%order%level_end)
    call this%levelsort(this%order, istat, root)
    if (istat /= 0_i8) call error_message("river%order: found cycle")
  end subroutine river_order

  !> \brief Calculate flow accumulation
  subroutine river_facc(this)
    use mo_message, only: error_message
    class(river_t), intent(inout) :: this
    integer(i8) :: i, j, n, m, k
    logical :: reverse_order
    ! integer(i8), pointer :: sources(:)
    if (allocated(this%facc)) deallocate(this%facc)
    if (.not.allocated(this%order%id)) call error_message("river%calc_facc: order not initialized")
    reverse_order = .not.this%order%to_root
    if (reverse_order) call this%order%reverse()
    allocate(this%facc(this%n_nodes))
    do i = 1_i8, this%order%n_levels
      !$omp parallel do default(shared) private(n,m,k) schedule(static)
      do j = this%order%level_start(i), this%order%level_end(i)
        n = this%order%id(j)
        this%facc(n) = 1_i4
        k = this%off_up(n)
        ! call this%src_view(n, sources)
        do m = k, k + this%n_up(n) - 1_i8
          this%facc(n) = this%facc(n) + this%facc(this%up(m))
        end do
      end do
      !$omp end parallel do
    end do
    ! restore original order if it was reversed
    if (reverse_order) call this%order%reverse()
  end subroutine river_facc

  !> \brief Calculate flow accumulation
  subroutine river_label_subcatchments(this, label_map, selected_nodes, labels, default_label)
    use mo_message, only: error_message
    class(river_t), intent(inout) :: this
    integer(i4), allocatable, intent(out) :: label_map(:) !< subcatchment labels
    integer(i8), dimension(:), intent(in) :: selected_nodes !< nodes to label subcatchments from
    integer(i4), dimension(:), optional, intent(in) :: labels !< labels used for subcatchments
    integer(i4), optional, intent(in) :: default_label !< default label for unlabeled nodes (default: 0)
    integer(i8) :: i, j, n
    integer(i4) :: def, m, k
    logical :: reverse_order

    if (present(labels)) then
      if (size(labels, kind=i8) /= size(selected_nodes, kind=i8)) then
        call error_message("river_label_subcatchments: size of labels does not match size of selected_nodes")
      end if
    end if

    if (.not.allocated(this%order%id)) call error_message("river%calc_facc: order not initialized")

    ! we need to traverse from downstream to upstream, so reverse order if it is currently from headwater to root
    reverse_order = this%order%to_root
    if (reverse_order) call this%order%reverse()

    def = optval(default_label, 0_i4)
    allocate(label_map(this%n_nodes))

    !$omp parallel do default(shared) schedule(static)
    do i = 1_i8, this%n_nodes
      label_map(i) = def
    end do
    !$omp end parallel do

    if (present(labels)) then
      do m = 1_i4, size(selected_nodes, kind=i4)
        if (labels(m) == def) then
          call error_message("river_label_subcatchments: given label equals default label")
        end if
        do k = m + 1_i4, size(selected_nodes, kind=i4)
          if (labels(m) == labels(k)) then
            call error_message("river_label_subcatchments: given label occurs more than once")
          end if
        end do
        if (label_map(selected_nodes(m)) /= def) then
          call error_message("river_label_subcatchments: selected node occurs more than once")
        end if
        label_map(selected_nodes(m)) = labels(m)
      end do
    else
      do m = 1_i4, size(selected_nodes, kind=i4)
        if (label_map(selected_nodes(m)) /= def) then
          call error_message("river_label_subcatchments: selected node occurs more than once")
        end if
        label_map(selected_nodes(m)) = m
      end do
    end if

    do i = 1_i8, this%order%n_levels
      !$omp parallel do default(shared) private(n) schedule(static)
      do j = this%order%level_start(i), this%order%level_end(i)
        n = this%order%id(j)
        if (this%down(n) == 0_i8) cycle ! sink
        if (label_map(n) /= def) cycle ! pre-labeled
        label_map(n) = label_map(this%down(n))
      end do
      !$omp end parallel do
    end do

    ! restore original order if it was reversed
    if (reverse_order) call this%order%reverse()

  end subroutine river_label_subcatchments

  !> \brief Calculate upstream area for each node (inclusive).
  subroutine river_upstream_area(this)
    use mo_message, only: error_message
    class(river_t), intent(inout), target :: this
    integer(i8) :: i, j, n, m, k
    logical :: reverse_order
    if (allocated(this%upstream_area)) deallocate(this%upstream_area)
    if (.not.allocated(this%order%id)) call error_message("river%calc_upstream_area: order not initialized")
    reverse_order = .not.this%order%to_root
    if (reverse_order) call this%order%reverse()
    allocate(this%upstream_area(this%n_nodes))
    do i = 1_i8, size(this%order%level_start, kind=i8)
      !$omp parallel do default(shared) private(n,m,k) schedule(static)
      do j = this%order%level_start(i), this%order%level_end(i)
        n = this%order%id(j)
        ! area fraction not allocated for d8 river
        if (this%scc) then
          this%upstream_area(n) = this%grid%cell_area(this%node_cell(n)) * this%area_fraction(n)
        else
          this%upstream_area(n) = this%grid%cell_area(n)
        end if
        k = this%off_up(n)
        do m = k, k + this%n_up(n) - 1_i8
          this%upstream_area(n) = this%upstream_area(n) + this%upstream_area(this%up(m))
        end do
      end do
      !$omp end parallel do
    end do
    ! restore original order if it was reversed
    if (reverse_order) call this%order%reverse()
  end subroutine river_upstream_area

  !> \brief Calculate D8 flow direction from network
  subroutine river_fdir(this)
    use mo_message, only: error_message
    class(river_t), intent(inout) :: this
    integer(i4) :: from(2), to(2)
    integer(i8) :: i, j
    integer(i4) :: dy ! direction of north in the grid matrix (1/-1)
    logical :: periodic ! periodic latlon grid
    if (this%scc) call error_message("river%calc_fdir: can not calculated fdir for a SCC river.")
    if (allocated(this%fdir)) deallocate(this%fdir)
    allocate(this%fdir(this%n_nodes))
    periodic = this%grid%is_periodic()
    dy = -1_i4 ! top-down grid starts north
    if (this%grid%y_direction==bottom_up) dy = 1_i4
    !$omp parallel do default(shared) private(j, from, to) schedule(static)
    do i = 1_i8, this%n_nodes
      j = this%down(i)
      if (j==0_i8) then
        this%fdir(i) = sink
        cycle ! sink
      end if
      from = this%grid%cell_ij(i,:)
      to = this%grid%cell_ij(j,:)
      ! the pathological case here is a periodic grid with 2 cells along lon-axis (*lol*)
      if (periodic .and. abs(from(1)-to(1))>1_i4) then
        ! set the node at the grid end to the begining
        if (from(1)==this%grid%nx) from(1) = 0_i4
        if (to(1)==this%grid%nx) to(1) = 0_i4
      end if
      this%fdir(i) = get_fdir(from, to, dy)
    end do
    !$omp end parallel do
  end subroutine river_fdir

  !> \brief Calculate link lengths
  subroutine river_length(this)
    use mo_constants, only: SQRT2_dp
    class(river_t), intent(inout) :: this
    integer(i8) :: i, to
    real(dp), allocatable, dimension(:) :: x_axis, y_axis
    if (this%scc) call error_message("river%calc_length: can't calculate D8 link lengths for a SCC river.")
    if (allocated(this%link_length)) deallocate(this%link_length)
    allocate(this%link_length(this%grid%ncells))
    if (this%grid%coordsys == cartesian) then
      !$omp parallel do default(shared) schedule(static)
      do i = 1_i8, this%grid%ncells
        if (this%down(i) == 0_i8) then
          ! sinks ain't links
          this%link_length(i) = 0.0_dp
          cycle
        end if
        select case (this%fdir(i))
          case(d8_E, d8_S, d8_W, d8_N)
            this%link_length(i) = this%grid%cellsize
          case(d8_SE, d8_SW, d8_NW, d8_NE)
            this%link_length(i) = SQRT2_dp * this%grid%cellsize
        end select
      end do
      !$omp end parallel do
    else
      x_axis = this%grid%x_axis()
      y_axis = this%grid%y_axis()
      !$omp parallel do default(shared) private(to) schedule(static)
      do i = 1_i8, this%grid%ncells
        to = this%down(i)
        if (to == 0_i8) then
          ! sinks ain't links
          this%link_length(i) = 0.0_dp
          cycle
        end if
        this%link_length(i) = dist_latlon( &
          lat1=y_axis(this%grid%cell_ij(i, 2)), &
          lon1=x_axis(this%grid%cell_ij(i, 1)), &
          lat2=y_axis(this%grid%cell_ij(to, 2)), &
          lon2=x_axis(this%grid%cell_ij(to, 1)))
      end do
      !$omp end parallel do
    end if
  end subroutine river_length

  !> \brief Calculate link slopes
  subroutine river_slope(this, dem)
    class(river_t), intent(inout) :: this
    real(dp), dimension(this%grid%ncells), intent(in) :: dem !< elevation [m]
    integer(i8) :: i, j, k
    if (this%scc) call error_message("river%calc_slope: can't calculate slope from DEM for a SCC river.")
    if (.not.allocated(this%link_length)) call error_message("river%calc_slope: can't calculate slope without link length.")
    if (allocated(this%link_slope)) deallocate(this%link_slope)
    allocate(this%link_slope(this%grid%ncells), source=0.0_dp)
    k = 0_i8
    !$omp parallel do default(shared) private(i,j)
    do i = 1_i8, this%grid%ncells
      j = this%down(i)
      if (j == 0_i8) cycle ! sinks ain't links
      if (dem(i) < dem(j)) then
        k = k + 1_i8
        ! *, "link is going upwards: ", i, j, dem(i), dem(j)
      end if
      this%link_slope(i) = max(( dem(i) - dem(j) ) / this%link_length(i), 0.0_dp) * 100.0_dp ! in %
    end do
    !$omp end parallel do
  end subroutine river_slope

  !> \brief Calculate the celerity
  !> \author Stephan Thober
  !> \author Matthias Kelbling
  !> \author Sebastian Müller
  subroutine river_celerity(this, gamma, constant_celerity, slope, mask)
    implicit none
    class(river_t), intent(inout) :: this
    real(dp), intent(in) :: gamma !< model parameter: c_i = gamma * sqrt(s_i) or c = gamma
    logical, optional, intent(in) :: constant_celerity !< whether celerity is assumed constant: c = gamma (default: .false.)
    real(dp), optional, intent(in) :: slope(:) !< [%] river slope, will be stored in link_slope if provided
    logical, optional, intent(in) :: mask(:) !< mask for slope smoothing (may come from upscaled river network)
    real(dp), allocatable :: smooth_slope(:)
    logical :: smoothing
    integer(i8) :: i

    ! constant celerity
    if (optval(constant_celerity, .false.)) then
      allocate(this%celerity(this%n_nodes))
      !$omp parallel do default(shared)
      do i = 1_i8, this%n_nodes
        this%celerity(i) = gamma
      end do
      !$omp end parallel do
      return
    end if

    ! need slope for non constant celerity
    if (present(slope)) then
      if (.not.allocated(this%link_slope)) allocate(this%link_slope(this%n_nodes))
      !$omp parallel do default(shared)
      do i = 1_i8, this%n_nodes
        this%link_slope(i) = slope(i)
      end do
      !$omp end parallel do
    end if

    if (.not.allocated(this%link_slope)) call error_message("river%calc_celerity: need slope to calculate celerity.")
    allocate(smooth_slope(this%n_nodes))
    !$omp parallel do default(shared)
    do i = 1_i8, this%n_nodes
      if (this%link_slope(i) < 0.1_dp) then
        ! set min val for river slope to enable h-mean
        smooth_slope(i) = 0.1_dp
      else
        smooth_slope(i) = this%link_slope(i)
      end if
    end do
    !$omp end parallel do

    ! smooth river slope if there is more than one cell
    smoothing = this%n_nodes > 1_i8
    if (present(mask)) smoothing = count(mask, kind=i8) > 1_i8
    if (smoothing) call river_smooth_slope(smooth_slope, mask)
    ! calculate celerity
    allocate(this%celerity(this%n_nodes))
    !$omp parallel do default(shared)
    do i = 1_i8, this%n_nodes
      this%celerity(i) = gamma * sqrt(smooth_slope(i) / 100.0_dp)
    end do
    !$omp end parallel do

  end subroutine river_celerity

  !> \brief Apply the hard-coded upper MAD slope filter used for river celerity.
  !> \details This replaces the generic MAD routine because its median selection uses default-kind indices.
  !! The implementation is equivalent to mad(arr, z=2.25_dp, mask=mask, tout="u", mval=0.1_dp).
  subroutine river_smooth_slope(slope, mask)
    real(dp), intent(inout) :: slope(:)
    logical, optional, intent(in) :: mask(:)
    integer(i8), parameter :: block_size = 1048576_i8
    integer(i8) :: n, nblocks, block, first, last, i, j, n_selected
    integer(i8), allocatable :: block_count(:), block_offset(:)
    real(dp), allocatable :: values(:)
    real(dp) :: median_slope, median_deviation, threshold
    logical :: fallback_to_floor

    n = size(slope, kind=i8)
    if (present(mask)) then
      if (size(mask, kind=i8) /= n) call error_message("river_smooth_slope: slope and mask have different size")
    end if
    nblocks = (n - 1_i8) / block_size + 1_i8
    allocate(block_count(nblocks), block_offset(nblocks))

    ! The generic MAD implementation excludes the floor value. If no values remain,
    ! it falls back to all floor values, independently of the input mask.
    !$omp parallel do default(shared) private(first, last, i) schedule(static)
    do block = 1_i8, nblocks
      first = (block - 1_i8) * block_size + 1_i8
      last = min(block * block_size, n)
      block_count(block) = 0_i8
      do i = first, last
        if (river_slope_selected(slope(i), i, mask, .false.)) block_count(block) = block_count(block) + 1_i8
      end do
    end do
    !$omp end parallel do
    fallback_to_floor = sum(block_count) == 0_i8

    !$omp parallel do default(shared) private(first, last, i) schedule(static)
    do block = 1_i8, nblocks
      first = (block - 1_i8) * block_size + 1_i8
      last = min(block * block_size, n)
      block_count(block) = 0_i8
      do i = first, last
        if (river_slope_selected(slope(i), i, mask, fallback_to_floor)) then
          block_count(block) = block_count(block) + 1_i8
        end if
      end do
    end do
    !$omp end parallel do

    n_selected = 0_i8
    do block = 1_i8, nblocks
      block_offset(block) = n_selected
      n_selected = n_selected + block_count(block)
    end do
    if (n_selected < 2_i8) call error_message("river_smooth_slope: need at least two slopes for MAD smoothing")
    allocate(values(n_selected))

    !$omp parallel do default(shared) private(first, last, i, j) schedule(static)
    do block = 1_i8, nblocks
      first = (block - 1_i8) * block_size + 1_i8
      last = min(block * block_size, n)
      j = block_offset(block)
      do i = first, last
        if (river_slope_selected(slope(i), i, mask, fallback_to_floor)) then
          j = j + 1_i8
          values(j) = slope(i)
        end if
      end do
    end do
    !$omp end parallel do

    median_slope = river_median_i8(values)
    !$omp parallel do default(shared)
    do i = 1_i8, n_selected
      values(i) = abs(values(i) - median_slope)
    end do
    !$omp end parallel do
    median_deviation = river_median_i8(values)
    threshold = median_deviation * 2.25_dp / 0.6745_dp
    deallocate(values, block_count, block_offset)

    !$omp parallel do default(shared)
    do i = 1_i8, n
      if (river_slope_selected(slope(i), i, mask, fallback_to_floor)) then
        if (slope(i) > median_slope + threshold) slope(i) = median_slope + threshold
      end if
    end do
    !$omp end parallel do
  end subroutine river_smooth_slope

  !> \brief Select the elements participating in the hard-coded river slope MAD filter.
  pure logical function river_slope_selected(value, index, mask, fallback_to_floor)
    real(dp), intent(in) :: value
    integer(i8), intent(in) :: index
    logical, optional, intent(in) :: mask(:)
    logical, intent(in) :: fallback_to_floor

    if (fallback_to_floor) then
      river_slope_selected = abs(value - 0.1_dp) < tiny(1.0_dp)
    else
      river_slope_selected = .not.(abs(value - 0.1_dp) < tiny(1.0_dp))
      if (present(mask)) river_slope_selected = river_slope_selected .and. mask(index)
    end if
  end function river_slope_selected

  !> \brief Return the median using 64-bit indices, rearranging values in place.
  function river_median_i8(values) result(median)
    real(dp), intent(inout) :: values(:)
    real(dp) :: median, previous
    integer(i8) :: n, middle

    n = size(values, kind=i8)
    if (mod(n, 2_i8) == 0_i8) then
      middle = n / 2_i8 + 1_i8
      median = river_nth_i8(values, middle)
      previous = maxval(values(:middle - 1_i8))
      median = 0.5_dp * (median + previous)
    else
      median = river_nth_i8(values, (n + 1_i8) / 2_i8)
    end if
  end function river_median_i8

  !> \brief Return the n-th smallest element using 64-bit indices, rearranging values in place.
  function river_nth_i8(values, n) result(nth)
    real(dp), intent(inout) :: values(:)
    integer(i8), intent(in) :: n
    real(dp) :: nth, pivot, tmp
    integer(i8) :: left, right, i, j

    left = 1_i8
    right = size(values, kind=i8)
    do while (left < right)
      pivot = values(n)
      i = left
      j = right
      do
        do while (values(i) < pivot)
          i = i + 1_i8
        end do
        do while (pivot < values(j))
          j = j - 1_i8
        end do
        if (i <= j) then
          tmp = values(i)
          values(i) = values(j)
          values(j) = tmp
          i = i + 1_i8
          j = j - 1_i8
        end if
        if (i > j) exit
      end do
      if (j < n) left = i
      if (n < i) right = j
    end do
    nth = values(n)
  end function river_nth_i8

  !> \brief Select values from sub-nodes for each cell from an array of values on nodes.
  function river_select_cell_values(this, values) result(select)
    class(river_t), intent(in) :: this
    real(dp), dimension(this%n_nodes), intent(in) :: values
    real(dp), dimension(this%grid%ncells):: select
    integer(i8) :: i
    !$omp parallel do default(shared)
    do i = 1_i8, this%grid%ncells
      select(i) = values(this%cell_node_select(i))
    end do
    !$omp end parallel do
  end function river_select_cell_values

  !> \brief Export river arrays to netcdf
  subroutine river_export(this, path, sub_map, leaving, stream_mask, stream_sub, highlight, factor)
    class(river_t), intent(in) :: this
    character(*), intent(in) :: path !< path to the file
    integer(i4), intent(in), optional :: sub_map(this%n_nodes) !< map of sub-catchment IDs
    logical, intent(in), optional :: leaving(this%n_nodes) !< mask of leaving cells
    logical, intent(in), optional :: stream_mask(this%n_nodes) !< stream mask
    integer(i4), intent(in), optional :: stream_sub(this%n_nodes) !< steam sub catchment ID
    logical, intent(in), optional :: highlight(:) !< highlight cells in leaving/stream_mask/stream_sub
    integer(i4), intent(in), optional :: factor !< upscaling factor for checkerboard pattern -1/0
    integer(i4), allocatable :: tmp(:)
    type(output_dataset) :: ds
    type(var), allocatable :: vars(:)
    integer(i8) :: i
    integer(i4) :: j, px, py
    integer(i4), allocatable :: level(:), pattern(:)
    if (this%scc) call error_message("river%export: can not export SCC river to gridded netcdf")
    if (present(factor)) then
      allocate(pattern(this%grid%ncells), source=0_i4)
      do i=1_i8, this%grid%ncells
        px = mod((this%grid%cell_ij(i, 1) - 1_i4) / factor, 2_i4)
        if (this%grid%y_direction==bottom_up) then
          py = mod((this%grid%cell_ij(i, 2) - 1_i4) / factor, 2_i4)
        else
          py = mod((this%grid%ny - this%grid%cell_ij(i, 2)) / factor, 2_i4)
        end if
        if (px+py==1_i4) pattern(i) = -1_i4
      end do
    end if
    allocate(vars(0))
    vars = [vars, var("fdir", "flow direction", dtype="i16", kind="i2", static=.true.)]
    if (allocated(this%facc)) vars = [vars, var("facc", "flow accumulation", dtype="i32", static=.true.)]
    if (allocated(this%order%id)) vars = [vars, var("level", "order level", dtype="i32", static=.true.)]
    if (allocated(this%upstream_area)) vars = [vars, var("upstream_area", "upstream area", units="m2", dtype="f64", static=.true.)]
    if (allocated(this%link_length)) vars = [vars, var("length", "link length", dtype="f64", static=.true.)]
    if (allocated(this%link_slope)) vars = [vars, var("slope", "link slope", dtype="f64", static=.true.)]
    if (allocated(this%celerity)) vars = [vars, var("celerity", "celerity", dtype="f64", static=.true.)]
    if (present(sub_map)) vars = [vars, var("scc", "scc catchment id", dtype="i32", static=.true.)]
    if (present(leaving)) vars = [vars, var("leaving", "leaving", dtype="i32", static=.true.)]
    if (present(stream_mask)) vars = [vars, var("stream", "stream mask", dtype="i32", static=.true.)]
    if (present(stream_sub)) vars = [vars, var("stream_sub", "stream sub catchment ID", dtype="i32", static=.true.)]
    call ds%init(path, this%grid, vars)
    call ds%update("fdir", this%fdir)
    if (allocated(this%facc)) call ds%update("facc", this%facc)
    if (allocated(this%order%id)) then
      allocate(level(this%grid%ncells))
      do j = 1_i4, size(this%order%level_start)
        do i = this%order%level_start(j), this%order%level_end(j)
          level(this%order%id(i)) = j
        end do
      end do
      call ds%update("level", level)
    end if
    if (allocated(this%upstream_area)) call ds%update("upstream_area", this%upstream_area)
    if (allocated(this%link_length)) call ds%update("length", this%link_length)
    if (allocated(this%link_slope)) call ds%update("slope", this%link_slope)
    if (allocated(this%celerity)) call ds%update("celerity", this%celerity)
    if (present(sub_map)) call ds%update("scc", sub_map)
    if (present(leaving)) then
      allocate(tmp(this%n_nodes), source=0_i4)
      if (present(factor)) tmp = pattern
      where(leaving) tmp = 1_i4
      if (present(highlight)) where(highlight) tmp = 2_i4
      call ds%update("leaving", tmp)
      deallocate(tmp)
    end if
    if (present(stream_mask)) then
      allocate(tmp(this%n_nodes), source=0_i4)
      if (present(factor)) tmp = pattern
      where(stream_mask) tmp = 1_i4
      if (present(highlight)) where(highlight) tmp = 2_i4
      call ds%update("stream", tmp)
      deallocate(tmp)
    end if
    if (present(stream_sub)) then
      if (present(factor)) then
        tmp = pattern
        where(stream_sub>0) tmp = stream_sub
        if (present(highlight)) where(highlight) tmp = maxval(stream_sub) + 1_i4
        call ds%update("stream_sub", tmp)
      else
        call ds%update("stream_sub", stream_sub)
      end if
    end if
    call ds%write()
    call ds%close()
    deallocate(vars)
  end subroutine river_export

  !> \brief Write river network to restart file
  subroutine river_to_restart_file(this, path, append)
    class(river_t), intent(in) :: this
    character(*), intent(in) :: path !< NetCDF file path
    logical, optional, intent(in) :: append !< whether NetCDF file should be opened in append mode (default .false.)
    type(NcDataset) :: nc
    character(1) :: fmode
    fmode = "w"
    if ( optval(append, .false.) ) fmode = "a"
    nc = NcDataset(path, fmode)
    call this%to_restart_dataset(nc)
    call nc%close()
  end subroutine river_to_restart_file

  !> \brief Write river network to restart dataset
  subroutine river_to_restart_dataset(this, nc)
    class(river_t), intent(in) :: this
    type(NcDimension) :: two_dim, node_dim, link_dim, dims(2), order_dim, sink_dim, xdim, ydim
    type(NcVariable) :: node_x_var, node_y_var, nc_var
    type(NcDataset) :: nc
    integer(i8) :: i
    integer(i8), allocatable :: dummy2d(:, :)
    integer(i2), allocatable :: dummy2d_i2(:, :)

    if ( this%grid%coordsys == cartesian ) then
      xdim = nc%getDimension("x")
      ydim = nc%getDimension("y")
    else
      xdim = nc%getDimension("lon")
      ydim = nc%getDimension("lat")
    end if

    ! river network
    if (this%n_nodes > int(huge(1_i4),i8)) call error_message("river%write_restart: river network too large to write to restart.")
    node_dim = nc%setDimension("node", int(this%n_nodes, i4))  ! only works if network is not to huge for i4
    sink_dim = nc%setDimension("sink", size(this%sinks))
    link_dim = nc%setDimension("link", size(this%up))
    two_dim = nc%setDimension("Two", 2_i4)

    ! 1D network topology following UGRID conventions
    nc_var = nc%setVariable("river", "i32", dims(:0)) ! mesh variable as scalar integer
    call nc_var%setAttribute("cf_role", "mesh_topology")
    call nc_var%setAttribute("long_name", "river network definition")
    call nc_var%setAttribute("topology_dimension", 1_i4)  ! 0 - only nodes, 1 - with links
    call nc_var%setAttribute("edge_node_connectivity", "links")
    if ( allocated(this%node_x) .and. allocated(this%node_y) ) then
      call nc_var%setAttribute("node_coordinates", "river_node_x river_node_y")
    end if

    ! links (from-to flow)
    nc_var = nc%setVariable("links", "i64", [two_dim, link_dim])
    call nc_var%setAttribute("cf_role", "edge_node_connectivity")
    call nc_var%setAttribute("long_name", "river links definition")
    call nc_var%setAttribute("start_index", 1_i4)  ! fortran indices starting with 1
    allocate(dummy2d(2_i4, size(this%up)))
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i8, this%n_nodes
      if (this%n_up(i) == 0_i8) cycle
      dummy2d(1_i4, this%off_up(i):this%off_up(i)+this%n_up(i)-1_i8) = this%up(this%off_up(i):this%off_up(i)+this%n_up(i)-1_i8)
      dummy2d(2_i4, this%off_up(i):this%off_up(i)+this%n_up(i)-1_i8) = i
    end do
    !$omp end parallel do
    call nc_var%setData(dummy2d)
    deallocate(dummy2d)

    ! down
    call message("writing down to restart file")
    nc_var = nc%setVariable("down", "i64", [node_dim])
    call nc_var%setAttribute("long_name", "downstream cell")
    call nc_var%setFillValue(nodata_i8)
    call nc_var%setAttribute("missing_value", nodata_i8)
    call nc_var%setData(this%down)

    ! sinks
    call message("writing sinks to restart file")
    nc_var = nc%setVariable("sinks", "i64", [sink_dim])
    call nc_var%setAttribute("long_name", "sinks")
    call nc_var%setFillValue(nodata_i8)
    call nc_var%setAttribute("missing_value", nodata_i8)
    call nc_var%setData(this%sinks)

    ! n_up
    call message("writing n_up to restart file")
    nc_var = nc%setVariable("n_up", "i64", [node_dim])
    call nc_var%setAttribute("long_name", "number of upstream nodes")
    call nc_var%setFillValue(nodata_i8)
    call nc_var%setAttribute("missing_value", nodata_i8)
    call nc_var%setData(this%n_up)

    ! off_up
    call message("writing off_up to restart file")
    nc_var = nc%setVariable("off_up", "i64", [node_dim])
    call nc_var%setAttribute("long_name", "offset to upstream nodes")
    call nc_var%setFillValue(nodata_i8)
    call nc_var%setAttribute("missing_value", nodata_i8)
    call nc_var%setData(this%off_up)

    ! is_sink
    ! don't need to write is_sink - can be derived from down

    ! tags
    ! don't need to write tags - will be standard either way (1..n)

    if ( allocated(this%node_x) .and. allocated(this%node_y) ) then
      ! coordinates
      node_x_var = nc%setVariable("river_node_x", "f64", [node_dim])
      node_y_var = nc%setVariable("river_node_y", "f64", [node_dim])
      if (this%grid%coordsys==cartesian) then
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
      call node_x_var%setData(this%node_x)
      call node_y_var%setData(this%node_y)
    end if

    ! fdir
    if ( allocated(this%fdir) ) then
      call message("writing fdir to restart file")
      nc_var = nc%setVariable("fdir", "i16", [xdim, ydim], shuffle=.true.)
      call nc_var%setAttribute("long_name", "flow direction")
      call nc_var%setFillValue(nodata_i2)
      call nc_var%setAttribute("missing_value", nodata_i2)
      allocate(dummy2d_i2(this%grid%nx, this%grid%ny))
      call this%grid%unpack_into(this%fdir, dummy2d_i2)
      call nc_var%setData(dummy2d_i2)
      deallocate(dummy2d_i2)
    end if

    ! facc
    if ( allocated(this%facc) ) then
      call message("writing facc to restart file")
      nc_var = nc%setVariable("facc", "i32", [node_dim])
      call nc_var%setAttribute("long_name", "flow accumulation")
      call nc_var%setFillValue(nodata_i4)
      call nc_var%setAttribute("missing_value", nodata_i4)
      call nc_var%setData(this%facc)
    end if

    ! upstream_area
    if ( allocated(this%upstream_area) ) then
      call message("writing upstream_area to restart file")
      nc_var = nc%setVariable("upstream_area", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "upstream area")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%upstream_area)
    end if

    ! link_length
    if ( allocated(this%link_length) ) then
      call message("writing link_length to restart file")
      nc_var = nc%setVariable("link_length", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "link length")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%link_length)
    end if

    ! link_slope
    if ( allocated(this%link_slope) ) then
      call message("writing link_slope to restart file")
      nc_var = nc%setVariable("link_slope", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "average slope of link")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%link_slope)
    end if

    ! celerity
    if ( allocated(this%celerity) ) then
      call message("writing celerity to restart file")
      nc_var = nc%setVariable("celerity", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "streamflow celerity")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%celerity)
    end if

    ! scc state
    nc_var = nc%setVariable("scc", "i8", dims(:0)) ! scalar integer to indicate scc river
    call nc_var%setAttribute("long_name", "scc river flag")
    call nc_var%setData(merge(1_i1, 0_i1, this%scc))
    call nc_var%setAttribute("flag_values", [0_i1, 1_i1])
    call nc_var%setAttribute("flag_meanings", "not_scc scc")

    ! node_cell
    if ( allocated(this%node_cell) ) then
      call message("writing node_cell to restart file")
      nc_var = nc%setVariable("node_cell", "i64", [node_dim])
      call nc_var%setAttribute("long_name", "map node to grid cell")
      call nc_var%setFillValue(nodata_i8)
      call nc_var%setAttribute("missing_value", nodata_i8)
      call nc_var%setData(this%node_cell)
    end if

    ! cell_node_select
    if ( allocated(this%cell_node_select) ) then
      call message("writing cell_node_select to restart file")
      nc_var = nc%setVariable("cell_node_select", "i64", [xdim, ydim], shuffle=.true.)
      call nc_var%setAttribute("long_name", "cell node select")
      call nc_var%setFillValue(nodata_i8)
      call nc_var%setAttribute("missing_value", nodata_i8)
      allocate(dummy2d(this%grid%nx, this%grid%ny))
      call this%grid%unpack_into(this%cell_node_select, dummy2d)
      call nc_var%setData(dummy2d)
      deallocate(dummy2d)
    end if

    ! area_fraction
    if ( allocated(this%area_fraction) ) then
      call message("writing area_fraction to restart file")
      nc_var = nc%setVariable("area_fraction", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "area fraction for each node in cell")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(this%area_fraction)
    end if

    ! order
    if ( allocated(this%order%id) ) then
      order_dim = nc%setDimension("order_dim", int(this%order%n_levels, i4))
      ! order%id
      nc_var = nc%setVariable("order_id", "i64", [node_dim])
      call nc_var%setAttribute("long_name", "id in order")
      call nc_var%setFillValue(nodata_i8)
      call nc_var%setAttribute("missing_value", nodata_i8)
      call nc_var%setData(this%order%id)
      ! order%level_start
      nc_var = nc%setVariable("order_level_start", "i64", [order_dim])
      call nc_var%setAttribute("long_name", "level start in order")
      call nc_var%setData(this%order%level_start)
      ! order%level_end
      nc_var = nc%setVariable("order_level_end", "i64", [order_dim])
      call nc_var%setAttribute("long_name", "level end in order")
      call nc_var%setData(this%order%level_end)
      ! order%level_size
      nc_var = nc%setVariable("order_level_size", "i64", [order_dim])
      call nc_var%setAttribute("long_name", "level size in order")
      call nc_var%setData(this%order%level_size)
    end if

  end subroutine river_to_restart_dataset

  !> \brief Read river network from restart file
  subroutine river_from_restart_file(this, path, grid)
    class(river_t), intent(inout) :: this
    character(*), intent(in) :: path !< NetCDF file path
    type(grid_t), pointer, intent(in) :: grid !< grid the river network is defined on
    type(NcDataset) :: nc
    nc = NcDataset(path, "r")
    call this%from_restart_dataset(nc, grid)
    call nc%close()
  end subroutine river_from_restart_file

  !> \brief Read river network from restart dataset
  subroutine river_from_restart_dataset(this, nc, grid)
    class(river_t), intent(inout) :: this
    type(NcDataset), intent(in) :: nc !< netcdf dataset to read from
    type(grid_t), pointer, intent(in) :: grid !< grid the river network is defined on
    type(NcVariable) :: nc_var
    type(NcDimension) :: nc_dim

    integer(i8) :: i, n_sinks, n_links
    integer(i1), allocatable :: dummy_i1
    integer(i2), allocatable :: dummy2di2(:, :)
    integer(i8), allocatable :: dummy2di8(:, :)

    ! reset all attributes
    call this%clean()

    ! grid should be read separately and passed as argument to avoid circular dependency between river and grid
    this%grid => grid

    nc_dim = nc%getDimension("node")
    this%n_nodes = int(nc_dim%getLength(), i8)
    nc_dim = nc%getDimension("sink")
    n_sinks = int(nc_dim%getLength(), i8)
    nc_dim = nc%getDimension("link")
    n_links = int(nc_dim%getLength(), i8)

    nc_var = nc%getVariable("links")
    allocate(this%up(n_links))
    call nc_var%getData(dummy2di8)
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i8, n_links
      this%up(i) = dummy2di8(1_i8, i)
    end do
    !$omp end parallel do
    deallocate(dummy2di8)

    nc_var = nc%getVariable("down")
    allocate(this%down(this%n_nodes))
    call nc_var%readInto(this%down)

    allocate(this%is_sink(this%n_nodes))
    allocate(this%tags(this%n_nodes))
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i8, this%n_nodes
      this%is_sink(i) = this%down(i) == 0_i8
      this%tags(i) = i
    end do
    !$omp end parallel do

    nc_var = nc%getVariable("sinks")
    allocate(this%sinks(n_sinks))
    call nc_var%readInto(this%sinks)

    nc_var = nc%getVariable("n_up")
    allocate(this%n_up(this%n_nodes))
    call nc_var%readInto(this%n_up)

    nc_var = nc%getVariable("off_up")
    allocate(this%off_up(this%n_nodes))
    call nc_var%readInto(this%off_up)

    if (nc%hasVariable("river_node_x")) then
      nc_var = nc%getVariable("river_node_x")
      allocate(this%node_x(this%n_nodes))
      call nc_var%readInto(this%node_x)
    end if

    if (nc%hasVariable("river_node_y")) then
      nc_var = nc%getVariable("river_node_y")
      allocate(this%node_y(this%n_nodes))
      call nc_var%readInto(this%node_y)
    end if

    if (nc%hasVariable("fdir")) then
      nc_var = nc%getVariable("fdir")
      call nc_var%getData(dummy2di2)
      allocate(this%fdir(this%grid%ncells))
      call this%grid%pack_into(dummy2di2, this%fdir)
      deallocate(dummy2di2)
    end if

    if (nc%hasVariable("facc")) then
      nc_var = nc%getVariable("facc")
      allocate(this%facc(this%n_nodes))
      call nc_var%readInto(this%facc)
    end if

    if (nc%hasVariable("upstream_area")) then
      nc_var = nc%getVariable("upstream_area")
      allocate(this%upstream_area(this%n_nodes))
      call nc_var%readInto(this%upstream_area)
    end if

    if (nc%hasVariable("link_length")) then
      nc_var = nc%getVariable("link_length")
      allocate(this%link_length(this%n_nodes))
      call nc_var%readInto(this%link_length)
    end if

    if (nc%hasVariable("link_slope")) then
      nc_var = nc%getVariable("link_slope")
      allocate(this%link_slope(this%n_nodes))
      call nc_var%readInto(this%link_slope)
    end if

    if (nc%hasVariable("celerity")) then
      nc_var = nc%getVariable("celerity")
      allocate(this%celerity(this%n_nodes))
      call nc_var%readInto(this%celerity)
    end if

    if (nc%hasVariable("node_cell")) then
      nc_var = nc%getVariable("node_cell")
      allocate(this%node_cell(this%n_nodes))
      call nc_var%readInto(this%node_cell)
    end if

    if (nc%hasVariable("cell_node_select")) then
      nc_var = nc%getVariable("cell_node_select")
      call nc_var%getData(dummy2di8)
      allocate(this%cell_node_select(this%grid%ncells))
      call this%grid%pack_into(dummy2di8, this%cell_node_select)
      deallocate(dummy2di8)
    end if

    if (nc%hasVariable("area_fraction")) then
      nc_var = nc%getVariable("area_fraction")
      allocate(this%area_fraction(this%n_nodes))
      call nc_var%readInto(this%area_fraction)
    end if

    if (nc%hasDimension("order_dim")) then
      this%order%to_root = .false.
      nc_dim = nc%getDimension("order_dim")
      this%order%n_levels = int(nc_dim%getLength(), i8)

      nc_var = nc%getVariable("order_id")
      allocate(this%order%id(this%n_nodes))
      call nc_var%readInto(this%order%id)

      nc_var = nc%getVariable("order_level_start")
      allocate(this%order%level_start(this%order%n_levels))
      call nc_var%readInto(this%order%level_start)

      nc_var = nc%getVariable("order_level_end")
      allocate(this%order%level_end(this%order%n_levels))
      call nc_var%readInto(this%order%level_end)

      nc_var = nc%getVariable("order_level_size")
      allocate(this%order%level_size(this%order%n_levels))
      call nc_var%readInto(this%order%level_size)
    end if

    this%scc = .false.
    if (nc%hasVariable("scc")) then
      nc_var = nc%getVariable("scc")
      call nc_var%getData(dummy_i1)
      this%scc = dummy_i1 == 1_i1
    end if

  end subroutine river_from_restart_dataset

  subroutine river_destroy(this)
    class(river_t), intent(inout) :: this
    if (allocated(this%fdir)) deallocate(this%fdir)
    if (allocated(this%facc)) deallocate(this%facc)
    if (allocated(this%is_sink)) deallocate(this%is_sink)
    if (allocated(this%upstream_area)) deallocate(this%upstream_area)
    if (allocated(this%link_length)) deallocate(this%link_length)
    if (allocated(this%link_slope)) deallocate(this%link_slope)
    if (allocated(this%celerity)) deallocate(this%celerity)
    if (allocated(this%node_x)) deallocate(this%node_x)
    if (allocated(this%node_y)) deallocate(this%node_y)
    if (allocated(this%node_cell)) deallocate(this%node_cell)
    if (allocated(this%cell_node_select)) deallocate(this%cell_node_select)
    if (allocated(this%area_fraction)) deallocate(this%area_fraction)
    call this%destroy()
  end subroutine river_destroy

end module mo_river
