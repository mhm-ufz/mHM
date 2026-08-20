!> \file    mo_river_upscaler.f90
!> \copydoc mo_river_upscaler

!> \brief   River upscaler.
!> \details This module contains an upscaler for river networks.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    May 2025
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_upscaler

  use mo_kind, only: i4, i8, dp
  use mo_constants, only: nodata_i4
  use mo_river, only: river_t, d8_leave_E, d8_leave_S, d8_leave_W, d8_leave_N, &
    d8_leave_SE, d8_leave_SW, d8_leave_NE, d8_leave_NW
  use mo_grid, only: grid_t, bottom_up, cartesian
  use mo_grid_scaler, only: scaler_t, down_scaling
  use mo_message, only: error_message, message
  use mo_utils, only: optval

  implicit none
  private

  integer(i4), public, parameter :: upscale_legacy = 0_i4 !< connect to the first entered coarse cell
  integer(i4), public, parameter :: upscale_flow = 1_i4 !< connect to the next downstream link endpoint

  !> \class river_upscaler_t
  !> \brief River network upscaler
  !> \details upscale river network respecting scc.
  !! Nodes are ordered in reverse sub-catchment order: first the nodes of the base catchment are added,
  !! then the nodes of the sub-catchment flowing into the base catchment and so on.
  type, public :: river_upscaler_t
    type(river_t), pointer :: fine_river => null() !< river definition at fine grid
    type(river_t), pointer :: coarse_river => null() !< river definition at coarse grid
    type(scaler_t) :: upscaler !< upscaler from fine to coarse grid
    ! coarse attributes derived from fine grid
    logical, allocatable :: leaving_cells(:) !< mask marking fine cells leaving a coarse cell size(fine\%ncells)
    logical, allocatable :: stream_mask(:) !< mask marking the upscaled stream at fine grid size(fine\%ncells)
    integer(i4), allocatable :: stream_sub(:) !< sub-catchment ID for stream cell at fine grid size(fine\%ncells)
    integer(i8), allocatable :: floodplain(:) !< map marking floodplain of a coarse river node with their id size(fine\%ncells)
    real(dp), allocatable :: floodplain_area(:) !< floodplain area for a coarse river node size(coarse\%n_nodes)
    logical, allocatable :: is_link_start(:) !< mask starting cell at fine grid of coarse river node size(fine\%n_nodes)
    integer(i8), allocatable :: link_start(:) !< starting cell at fine grid of coarse river node (0 for sink) size(coarse\%n_nodes)
    integer(i8), allocatable :: link_end(:) !< downstream fine endpoint of coarse river node (0 for sink) size(coarse\%n_nodes)
    integer(i8), allocatable :: sink_map(:) !< id of defining fine cell for a coarse sink (0 for no sink) size(coarse\%n_nodes)
    ! scc related attributes
    integer(i4), allocatable :: scc_map(:) !< sub catchment id map (id 'nsub' indicates base-catchment) size(fine\%n_nodes)
    integer(i8), allocatable :: scc_gauges(:) !< id of gauge for sub catchment on fine grid size(nsub-1)
    integer(i8), allocatable :: scc_coarse_gauges(:) !< id of gauge for sub catchment on coarse grid size(nsub-1)
    logical, allocatable :: is_scc_gauge(:) !< flag if node is a scc gauge size(fine%\n_nodes)
    logical, allocatable :: is_scc_coarse_gauge(:) !< flag if coarse node is a scc gauge size(coarse%\n_nodes)
    ! sub-catchment related attributes
    integer(i4) :: nsub = 1_i4 !< number of sub catchments
    logical, allocatable :: has_sub(:,:) !< flag if coarse cell contains node of a sub-catchment size(coarse\%ncells,nsub)
    integer(i4), allocatable :: n_sub_nodes(:) !< number of scc nodes per coarse grid cell size(coarse\%ncells)
    integer(i4), allocatable :: node_sub(:) !< map coarse node to sub catchment id size(coarse\%n_nodes)
    integer(i8), allocatable :: sub_size(:) !< number of coarse river nodes in each sub catchment size(nsub)
    integer(i8), allocatable :: sub_cum_size(:) !< cumulative sum of sub_size with (number of previous sub catchments nodes) size(nsub)
  contains
    procedure, public :: init => river_upscaler_init
    procedure, private :: init_scc_gauges => river_upscaler_init_scc_gauges
    procedure, private :: find_leaving_cells => river_upscaler_leaving
    procedure, private :: upscale => river_upscaler_upscale
    procedure, public :: calc_celerity => river_upscaler_celerity
    procedure, public :: node_from_cell_sub => river_upscaler_cell_sub
    procedure, public :: destroy => river_upscaler_destroy
  end type river_upscaler_t

contains

  !> \brief Setup river upscaler from fine river and coarse target grid.
  subroutine river_upscaler_init(this, &
    fine_river, coarse_river, coarse_grid, scc_gauges, scc_latlon, upscale_mode, length_percentile, tol)
    implicit none
    class(river_upscaler_t), intent(inout) :: this
    type(river_t), pointer, intent(in) :: fine_river !< river definition at fine grid
    type(river_t), pointer, intent(in) :: coarse_river !< pointer to coarse river definition to be determined
    type(grid_t), pointer, intent(in) :: coarse_grid !< coarse grid for the upscaled river network
    real(dp), dimension(:,:), intent(in), optional :: scc_gauges !< gauge locations on fine river dim 1: id, dim 2: (x,y)
    logical, optional, intent(in) :: scc_latlon !< Whether the scc gauges coordinates are geographical (default: .true.)
    integer(i4), optional, intent(in) :: upscale_mode !< upscaling mode (0: legacy, 1: FLOW-like; default: 1)
    real(dp), optional, intent(in) :: length_percentile !< [%] percentile for lower cut-off of upscaled link-length (40 by default)
    real(dp), optional, intent(in) :: tol !< tolerance for cell factor comparison (default: 1.e-7)
    integer(i4) :: mode
    real(dp) :: len_percentile

    mode = optval(upscale_mode, upscale_flow)
    len_percentile = optval(length_percentile, 40.0_dp)
    if (mode /= upscale_legacy .and. mode /= upscale_flow) then
      call error_message("river_upscaler: upscale_mode needs to be 0 (legacy) or 1 (FLOW-like)")
    end if
    if (len_percentile < 0.0_dp .or. len_percentile > 100.0_dp) then
      call error_message("river_upscaler: length_percentile needs to be in the range [0, 100]")
    end if

    this%fine_river => fine_river
    this%coarse_river => coarse_river
    this%coarse_river%grid => coarse_grid
    call this%upscaler%init(this%fine_river%grid, this%coarse_river%grid, tol=tol)
    if (this%upscaler%scaling_mode == down_scaling) call error_message("river_upscaler: target grid needs to be coarser then input")
    ! TODO: shortcut for same resolution of fine and coarse river

    ! sanity check to have facc available for upscaling
    if (.not.allocated(this%fine_river%facc)) then
      if (.not.allocated(this%fine_river%order%id)) call this%fine_river%calc_order(root=.true.)
      call this%fine_river%calc_facc()
    end if

    ! initialize scc related variables
    call message("river_upscaler: initialize scc")
    call this%init_scc_gauges(scc_gauges, scc_latlon)

    ! find leaving fine cells in every coarse cell
    call message("river_upscaler: find leaving cells")
    call this%find_leaving_cells()

    ! upscale graph and calculate stream features in one fine-river trace
    call message("river_upscaler: upscale river graph and stream features")
    call this%upscale(mode, len_percentile)

  end subroutine river_upscaler_init

  !> \brief Initialize SCC related variables
  subroutine river_upscaler_init_scc_gauges(this, scc_gauges, scc_latlon)
    class(river_upscaler_t), intent(inout) :: this
    real(dp), dimension(:,:), intent(in), optional :: scc_gauges !< gauge locations on fine river dim 1: id, dim 2: (x,y)
    logical, optional, intent(in) :: scc_latlon !< Whether the scc gauges coordinates are geographical (default: .true.)
    ! sub-catchment related attributes
    integer(i4) :: i, n
    integer(i8) :: k
    logical :: aux, latlon

    if (.not.present(scc_gauges)) then
      call message("river_upscaler%init_scc: no scc gauges provided, initialize without scc")
      this%nsub = 1_i4
      allocate(this%scc_gauges(0))
      allocate(this%scc_coarse_gauges(0))
      allocate(this%scc_map(this%fine_river%n_nodes))
      allocate(this%is_scc_gauge(this%fine_river%n_nodes))
      ! set default
      !$omp parallel do default(shared)
      do k = 1_i8, this%fine_river%n_nodes
        this%scc_map(k) = 1_i4
        this%is_scc_gauge(k) = .false.
      end do
      !$omp end parallel do
      this%coarse_river%scc = .false.
      return
    end if

    latlon = optval(scc_latlon, .true.)  ! assume geographical coordinates of stations by default
    if (latlon .and. (.not.this%fine_river%grid%has_aux_coords()) .and. this%fine_river%grid%coordsys == cartesian) then
      call error_message("river_upscaler%init_scc: gauge coordinates are geographical but grid has no lat-lon coordinates.")
    end if
    aux = this%fine_river%grid%has_aux_coords() .and. latlon

    n = size(scc_gauges, dim=1)
    this%nsub = n + 1_i4
    if (this%fine_river%scc) call error_message("river_upscaler%init_scc: need a D8 river to initialize SCC")
    if (.not.allocated(this%fine_river%facc)) call error_message("river_upscaler%init_scc: facc not available")

    ! find gauge cell ids
    allocate(this%scc_gauges(n))
    call message("river_upscaler%init_scc: find gauge cell ids")
    if (aux) then
      call message("river_upscaler%init_scc: find gauge cell with aux coordinates (slow)")
      !$omp parallel do default(shared)
      do i = 1_i4, n
        this%scc_gauges(i) = this%fine_river%grid%closest_cell_id(scc_gauges(i,:), use_aux=aux)
      end do
      !$omp end parallel do
    else
      this%scc_gauges = this%fine_river%grid%closest_cell_id_by_axes(scc_gauges)
    end if

    ! calculate scc map
    if (.not.allocated(this%fine_river%order%id)) call this%fine_river%calc_order(root=.true.)
    call this%fine_river%label_subcatchments(this%scc_map, this%scc_gauges, default_label=this%nsub) ! base catchment gets id n+1

    ! determine gauges
    allocate(this%is_scc_gauge(this%fine_river%n_nodes))
    call message("river_upscaler%init_scc: determine scc gauges on fine grid")
    ! set default
    !$omp parallel do default(shared)
    do k = 1_i8, this%fine_river%n_nodes
      this%is_scc_gauge(k) = .false.
    end do
    !$omp end parallel do
    ! fill in gauges
    !$omp parallel do default(shared)
    do i = 1_i4, n
      this%is_scc_gauge(this%scc_gauges(i)) = .true.
    end do
    !$omp end parallel do
    this%coarse_river%scc = this%nsub > 1_i4
  end subroutine river_upscaler_init_scc_gauges

  !> \brief Find all leaving cells of fine grid on a coarse cell border.
  subroutine river_upscaler_leaving(this)
    implicit none
    class(river_upscaler_t), intent(inout) :: this
    integer(i8), allocatable :: cells(:,:)
    integer(i4), allocatable :: scc_map(:,:)
    integer(i8) :: i
    integer(i4) :: j, ix, iy
    integer(i4) :: yl, yu, xl, xu ! lower and upper bounds for x and y
    integer(i4) :: yn, ys ! north and south y-bound depending on grid y-direction
    logical :: has_sub_init

    ! find all leaving cells (fine cells at coarse cell borders pointing out the coarse cell)
    !------------------------------------------------------------------
    !                            xl     xu
    !                        yn  NW--N--NE       sides:   E, S, W, N
    !                            |       |       corners: SE,SW,NW,NE
    !                            W  (i)  E       yu and yl with inverse sky direction if grid is top-down -> yn and ys
    !                            |       |
    !                        ys  SW--S--SE
    !------------------------------------------------------------------

    allocate(cells(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%gen_id_matrix(cells)

    allocate(this%leaving_cells(this%fine_river%grid%ncells))
    !$omp parallel do default(shared)
    do i = 1_i8, this%fine_river%grid%ncells
      this%leaving_cells(i) = .false.
    end do
    !$omp end parallel do

    allocate(this%has_sub(this%coarse_river%grid%ncells, this%nsub))
    has_sub_init = .not.this%coarse_river%scc
    !$omp parallel do default(shared)
    do j = 1_i4, this%nsub
      this%has_sub(:,j) = has_sub_init
    end do
    !$omp end parallel do

    if (this%coarse_river%scc) then
      allocate(scc_map(this%fine_river%grid%nx, this%fine_river%grid%ny))
      call this%fine_river%grid%unpack_into(this%scc_map, scc_map)
    end if

    !$omp parallel do default(shared) private(i, yn, ys, ix, iy, j, yl, yu, xl, xu)
    do i = 1_i8, this%coarse_river%grid%ncells
      call this%upscaler%coarse_bounds(i, xl, xu, yl, yu)
      ! determine north and south bound depending on y-direction
      if (this%fine_river%grid%y_direction == bottom_up) then
        yn = yu ! north y-bound is upper bound
        ys = yl ! south y-bound is lower bound
      else ! top_down
        yn = yl ! north y-bound is lower bound
        ys = yu ! south y-bound is upper bound
      end if

      ! searching on side E
      do iy = yl + 1_i4, yu - 1_i4
        if (.not.this%fine_river%grid%mask(xu,iy)) cycle
        if (any(this%fine_river%fdir(cells(xu,iy))==d8_leave_E)) this%leaving_cells(cells(xu,iy)) = .true.
      end do
      ! searching on side S
      do ix = xl + 1_i4, xu - 1_i4
        if (.not.this%fine_river%grid%mask(ix,ys)) cycle
        if (any(this%fine_river%fdir(cells(ix,ys))==d8_leave_S)) this%leaving_cells(cells(ix,ys)) = .true.
      end do
      ! searching on side W
      do iy = yl + 1_i4, yu - 1_i4
        if (.not.this%fine_river%grid%mask(xl,iy)) cycle
        if (any(this%fine_river%fdir(cells(xl,iy))==d8_leave_W)) this%leaving_cells(cells(xl,iy)) = .true.
      end do
      ! searching on side N
      do ix = xl + 1_i4, xu - 1_i4
        if (.not.this%fine_river%grid%mask(ix,yn)) cycle
        if (any(this%fine_river%fdir(cells(ix,yn))==d8_leave_N)) this%leaving_cells(cells(ix,yn)) = .true.
      end do
      ! searching on corner SE
      if (this%fine_river%grid%mask(xu,ys)) then
        if (any(this%fine_river%fdir(cells(xu,ys))==d8_leave_SE)) this%leaving_cells(cells(xu,ys)) = .true.
      end if
      ! searching on corner SW
      if (this%fine_river%grid%mask(xl,ys)) then
        if (any(this%fine_river%fdir(cells(xl,ys))==d8_leave_SW)) this%leaving_cells(cells(xl,ys)) = .true.
      end if
      ! searching on corner NE
      if (this%fine_river%grid%mask(xu,yn)) then
        if (any(this%fine_river%fdir(cells(xu,yn))==d8_leave_NE)) this%leaving_cells(cells(xu,yn)) = .true.
      end if
      ! searching on corner NW
      if (this%fine_river%grid%mask(xl,yn)) then
        if (any(this%fine_river%fdir(cells(xl,yn))==d8_leave_NW)) this%leaving_cells(cells(xl,yn)) = .true.
      end if

      ! count number of scc nodes
      if (this%coarse_river%scc) then
        do j = 1_i4, this%nsub
          if (any(scc_map(xl:xu,yl:yu)==j)) this%has_sub(i, j) = .true.
        end do
      end if
    end do
    !$omp end parallel do
  end subroutine river_upscaler_leaving

  !> \brief Setup the coarse graph and stream features with optional SCC nodes.
  !> \details Both modes trace source-inclusive and endpoint-exclusive fine links.
  subroutine river_upscaler_upscale(this, upscale_mode, length_percentile)
    use mo_percentile, only: percentile
    use mo_utils, only: prefix_sum
    implicit none
    class(river_upscaler_t), intent(inout) :: this
    integer(i4), intent(in) :: upscale_mode
    real(dp), intent(in) :: length_percentile
    integer(i8), allocatable :: down(:), ids(:,:)
    integer(i4), allocatable :: facc(:,:), scc_map(:,:)
    logical, allocatable :: leave_mask(:,:)
    integer(i8) :: i, k, node, next, cell, facc_max_i, n_nodes, n_links
    integer(i4) :: j, sub, facc_max, ix, iy, loc(2)
    integer(i4) :: yl, yu, xl, xu
    real(dp) :: length_cutoff

    if (.not.allocated(this%fine_river%facc)) call error_message("river_upscaler%upscale: facc not available")
    if (.not.allocated(this%fine_river%link_length)) call error_message("river_upscaler%upscale: link length not available")
    if (.not.allocated(this%fine_river%node_x)) call error_message("river_upscaler%upscale: node location not available")
    if (.not.allocated(this%fine_river%node_y)) call error_message("river_upscaler%upscale: node location not available")

    call message("river_upscaler%upscale: allocate coarse river attributes")
    allocate(this%sub_size(this%nsub))
    allocate(this%sub_cum_size(this%nsub))
    allocate(this%n_sub_nodes(this%coarse_river%grid%ncells))
    !$omp parallel do default(shared)
    do i = 1_i8, this%coarse_river%grid%ncells
      this%n_sub_nodes(i) = count(this%has_sub(i,:), kind=i4)
    end do
    !$omp end parallel do
    !$omp parallel do default(shared)
    do j = 1_i4, this%nsub
      this%sub_size(j) = count(this%has_sub(:,j), kind=i8)
    end do
    !$omp end parallel do
    call prefix_sum(this%sub_size, this%sub_cum_size, shift=1_i8, start=0_i8)
    n_nodes = this%sub_cum_size(this%nsub) + this%sub_size(this%nsub) ! total number of coarse nodes

    ! coarse river attributes
    if (this%coarse_river%scc) then
      allocate(this%coarse_river%node_cell(n_nodes))
      allocate(this%coarse_river%area_fraction(n_nodes))
    end if
    allocate(this%coarse_river%is_sink(n_nodes))
    allocate(this%coarse_river%link_length(n_nodes))
    allocate(this%coarse_river%node_x(n_nodes))
    allocate(this%coarse_river%node_y(n_nodes))
    ! upscaler attributes
    allocate(this%is_scc_coarse_gauge(n_nodes))
    allocate(this%node_sub(n_nodes))
    allocate(this%sink_map(n_nodes))
    allocate(this%link_start(n_nodes))
    allocate(this%link_end(n_nodes))

    ! initialize attributes
    !$omp parallel do default(shared)
    do k = 1_i8, n_nodes
      this%coarse_river%link_length(k) = 0.0_dp
      this%coarse_river%is_sink(k) = .false.
      this%is_scc_coarse_gauge(k) = .false.
      this%sink_map(k) = 0_i8
      this%link_start(k) = 0_i8
      this%link_end(k) = 0_i8
    end do
    !$omp end parallel do

    ! determine sub-catchment for each coarse node
    call message("river_upscaler%upscale: determine sub-catchment for each coarse node")
    !$omp parallel do default(shared) private(k, i)
    do j = 1_i4, this%nsub
      k = sum(this%sub_size(1_i4:j-1_i4))
      do i = 1_i8, this%coarse_river%grid%ncells
        if (.not.this%has_sub(i, j)) cycle
        k = k + 1_i8
        if (this%coarse_river%scc) this%coarse_river%node_cell(k) = i
        this%node_sub(k) = j
      end do
    end do
    !$omp end parallel do

    ! determine coarse gauge for each sub-catchment gauge
    call message("river_upscaler%upscale: determine coarse gauge for each sub-catchment gauge")
    allocate(this%scc_coarse_gauges(size(this%scc_gauges)))
    !$omp parallel do default(shared) private(k)
    do j = 1_i4, size(this%scc_gauges)
      k = this%node_from_cell_sub(this%upscaler%id_map(this%scc_gauges(j)), j)
      this%scc_coarse_gauges(j) = k
      this%is_scc_coarse_gauge(k) = .true.
    end do
    !$omp end parallel do

    ! unpack scc map for fine river
    allocate(scc_map(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%unpack_into(this%scc_map, scc_map)

    ! determine area fraction of each coarse node
    if (this%coarse_river%scc) then
      call message("river_upscaler%upscale: determine area fraction of each coarse node")
      !$omp parallel do default(shared)
      do k = 1_i8, n_nodes
        this%coarse_river%area_fraction(k) = this%upscaler%cell_fraction( &
          class_map   = scc_map, &
          coarse_cell = this%coarse_river%node_cell(k), &
          class_id    = this%node_sub(k))
      end do
      !$omp end parallel do
    end if

    ! unpack facc for fine river
    allocate(facc(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%unpack_into(this%fine_river%facc, facc)

    ! determine coarse sinks and their defining fine cell
    call message("river_upscaler%upscale: determine coarse sinks and their defining fine cell")
    !$omp parallel do default(shared) private(i, k, sub, yl, yu, xl, xu, node)
    do j = 1_i4, size(this%fine_river%sinks)
      i = this%fine_river%sinks(j)
      k = this%upscaler%id_map(i)
      sub = this%scc_map(i)
      ! sinks in base catchment are only sinks if they are the maximum facc in the coarse cell
      if (sub == this%nsub) then
        call this%upscaler%coarse_bounds(k, xl, xu, yl, yu)
        if (this%fine_river%facc(i) < maxval(facc(xl:xu,yl:yu), mask=(scc_map(xl:xu,yl:yu)==sub))) cycle
      end if
      node = this%node_from_cell_sub(k, sub)
      this%coarse_river%is_sink(node) = .true.
      this%sink_map(node) = i
    end do
    !$omp end parallel do

    ! unpack leaving cells for fine river
    allocate(leave_mask(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%unpack_into(this%leaving_cells, leave_mask)
    ! create ID matrix for fine river
    allocate(ids(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%gen_id_matrix(ids)

    ! determine starting fine cell for each coarse node
    call message("river_upscaler%upscale: determine starting fine cell for each coarse node")
    !$omp parallel do default(shared) private(k, sub, xl, xu, yl, yu, loc, ix, iy)
    do i = 1_i8, n_nodes
      if (this%coarse_river%is_sink(i)) cycle
      sub = this%node_sub(i)
      if (this%is_scc_coarse_gauge(i)) then
        this%link_start(i) = this%scc_gauges(sub)
      else
        if (this%coarse_river%scc) then
          k = this%coarse_river%node_cell(i)
          call this%upscaler%coarse_bounds(k, xl, xu, yl, yu)
          loc = maxloc(facc(xl:xu,yl:yu), mask=(leave_mask(xl:xu,yl:yu).and.(scc_map(xl:xu,yl:yu)==sub)))
        else
          k = i
          call this%upscaler%coarse_bounds(k, xl, xu, yl, yu)
          loc = maxloc(facc(xl:xu,yl:yu), mask=(leave_mask(xl:xu,yl:yu)))
        end if
        ix = xl + loc(1) - 1_i4
        iy = yl + loc(2) - 1_i4
        this%link_start(i) = ids(ix, iy)
      end if
    end do
    !$omp end parallel do
    deallocate(ids, facc, scc_map, leave_mask)

    ! initialize stream mask and stream sub-catchment map
    call message("river_upscaler%upscale: initialize stream mask and stream sub-catchment map")
    allocate(this%is_link_start(this%fine_river%n_nodes))
    allocate(this%stream_mask(this%fine_river%n_nodes))
    allocate(this%stream_sub(this%fine_river%n_nodes))
    !$omp parallel do default(shared)
    do i = 1_i8, this%fine_river%n_nodes
      this%is_link_start(i) = .false.
      this%stream_mask(i) = .false.
      this%stream_sub(i) = nodata_i4
    end do
    !$omp end parallel do

    ! initialize coarse node locations and mark link start cells
    call message("river_upscaler%upscale: initialize coarse node locations and mark link start cells")
    !$omp parallel do default(shared)
    do i = 1_i8, n_nodes
      if (this%coarse_river%is_sink(i)) then
        this%coarse_river%node_x(i) = this%fine_river%node_x(this%sink_map(i))
        this%coarse_river%node_y(i) = this%fine_river%node_y(this%sink_map(i))
      else
        this%coarse_river%node_x(i) = this%fine_river%node_x(this%link_start(i))
        this%coarse_river%node_y(i) = this%fine_river%node_y(this%link_start(i))
        this%is_link_start(this%link_start(i)) = .true.
      end if
    end do
    !$omp end parallel do

    ! mark stream mask and determine coarse link lengths and downstream nodes
    call message("river_upscaler%upscale: mark stream mask and determine coarse link lengths and downstream nodes")
    allocate(down(n_nodes))
    !$omp parallel do default(shared) private(cell, next, sub)
    do i = 1_i8, n_nodes
      if (this%coarse_river%is_sink(i)) then
        this%stream_mask(this%sink_map(i)) = .true.
        down(i) = 0_i8
        cycle
      end if
      cell = this%link_start(i)
      next = this%fine_river%down(cell)
      if (upscale_mode == upscale_legacy) then
        sub = this%scc_map(next)
        down(i) = this%node_from_cell_sub(this%upscaler%id_map(next), sub)
      end if
      do
        ! Concurrent traces may set the same fine cell to true; the result is invariant.
        this%stream_mask(cell) = .true.
        this%coarse_river%link_length(i) = this%coarse_river%link_length(i) + this%fine_river%link_length(cell)
        cell = this%fine_river%down(cell)
        ! if (this%leaving_cells(cell)) exit ! stop when leaving next coarse cell
        if (this%is_link_start(cell)) exit
        if (this%is_scc_gauge(cell)) exit
        if (this%fine_river%is_sink(cell)) exit
      end do
      this%stream_mask(cell) = .true.
      this%link_end(i) = cell
      if (upscale_mode == upscale_flow) then
        sub = this%scc_map(cell)
        down(i) = this%node_from_cell_sub(this%upscaler%id_map(cell), sub)
      end if
    end do
    !$omp end parallel do

    if (allocated(this%coarse_river%fdir)) deallocate(this%coarse_river%fdir)
    call message("river_upscaler: initialize coarse river branching DAG")
    call this%coarse_river%init(down)
    deallocate(down)

    if (this%coarse_river%scc) then
      call message("river_upscaler: find representative node for each coarse cell")
      allocate(this%coarse_river%cell_node_select(this%coarse_river%grid%ncells))
      !$omp parallel do default(shared) private(j, k, node, facc_max, facc_max_i)
      do i = 1_i8, this%coarse_river%grid%ncells
        facc_max = 0_i4
        facc_max_i = 0_i8
        do j = 1_i4, this%nsub
          k = this%node_from_cell_sub(i, j)
          if (k == 0_i8) cycle
          if (this%coarse_river%is_sink(k)) then
            node = this%sink_map(k)
          else
            node = this%link_start(k)
          end if
          if (this%fine_river%facc(node) > facc_max) then
            facc_max = this%fine_river%facc(node)
            facc_max_i = k
          end if
        end do
        this%coarse_river%cell_node_select(i) = facc_max_i
      end do
      !$omp end parallel do
    end if

    ! apply lower cut-off for coarse link lengths
    n_links = n_nodes - this%coarse_river%n_roots()
    if (n_links > 0_i8) then
      call message("river_upscaler: apply lower cut-off for coarse link lengths")
      length_cutoff = percentile(pack(this%coarse_river%link_length, mask=.not.this%coarse_river%is_sink), length_percentile)
      !$omp parallel do default(shared)
      do i = 1_i8, n_nodes
        if (this%coarse_river%is_sink(i)) cycle
        if (this%coarse_river%link_length(i) < length_cutoff) this%coarse_river%link_length(i) = length_cutoff
      end do
      !$omp end parallel do
    end if

    ! mark stream sub-catchment map for fine river
    call message("river_upscaler: mark stream sub-catchment map for fine river")
    !$omp parallel do default(shared)
    do i = 1_i8, this%fine_river%n_nodes
      if (this%stream_mask(i)) this%stream_sub(i) = this%scc_map(i)
    end do
    !$omp end parallel do

  end subroutine river_upscaler_upscale

  !> \brief calculate the celerity c_i from slope s_i (i - cell index)
  subroutine river_upscaler_celerity(this, gamma, constant_celerity, slope)
    implicit none
    class(river_upscaler_t), intent(inout) :: this
    real(dp), intent(in) :: gamma !< model parameter: c_i = gamma * sqrt(s_i) or c = gamma
    logical, optional, intent(in) :: constant_celerity !< whether celerity is assumed constant: c = gamma (default: .false.)
    real(dp), optional, intent(in) :: slope(:) !< [%] river slope on fine grid: size(fine\%ncells)
    integer(i8) :: i, cell
    real(dp) :: n
    logical :: constant

    ! first calculate celerity on fine river
    call this%fine_river%calc_celerity(gamma, constant_celerity, slope, this%stream_mask)

    constant = optval(constant_celerity, .false.)
    allocate(this%coarse_river%celerity(this%coarse_river%n_nodes))
    if (constant) then
      call message("river_upscaler: constant celerity assumed, set c_i = gamma for all coarse nodes")
      !$omp parallel do default(shared)
      do i = 1_i8, this%coarse_river%n_nodes
        this%coarse_river%celerity(i) = gamma
      end do
      !$omp end parallel do
      return
    end if

    call message("river_upscaler: calculate celerity on coarse river from fine river")
    !$omp parallel do default(shared) private(i, cell, n)
    do i = 1_i8, this%coarse_river%n_nodes
      if (this%coarse_river%is_sink(i)) then
        this%coarse_river%celerity(i) = 1.0_dp
        cycle
      end if
      cell = this%link_start(i)
      this%coarse_river%celerity(i) = 0.0_dp
      n = 0.0_dp
      ! one pass algorithm for harmonic mean:
      ! 0. M  = 0               -> 0 as initial value for the [M]ean of inverses
      ! 1. M += (1/v - M) / n   -> updating the mean with the weighted deviation of new value (n as counter)
      ! 2. H  = 1 / M           -> [H]armonic mean is then the inverse of M
      do while (cell /= this%link_end(i))
        n = n + 1.0_dp
        this%coarse_river%celerity(i) = this%coarse_river%celerity(i) &
          + ( 1.0_dp / this%fine_river%celerity(cell) - this%coarse_river%celerity(i) ) / n
        cell = this%fine_river%down(cell)
      end do
      ! finalize harmonic mean for celerity
      this%coarse_river%celerity(i) = 1.0_dp / this%coarse_river%celerity(i)
    end do
    !$omp end parallel do

  end subroutine river_upscaler_celerity

  !> \brief determine node ID from coarse cell and sub-catchment id
  !> \details Nodes are ordered by sub-catchments:
  !! 1. all coarse cells touching sub-catchment 1
  !! 2. all coarse cells touching sub-catchment 2 ...
  !!
  !! last: all cells from the base-catchment (not touched by any specified sub-catchment)
  pure function river_upscaler_cell_sub(this, cell, sub) result(id)
    implicit none
    class(river_upscaler_t), intent(in) :: this
    integer(i8), intent(in) :: cell !< cell id on coarse river grid (1..coarse\%ncells)
    integer(i4), intent(in) :: sub !< sub-catchment id from scc (1..nsub)
    integer(i8) :: id
    if (.not.this%has_sub(cell, sub)) then
      ! given cell doesn't contain the specified sub-catchment
      id = 0_i8
      return
    end if
    id = this%sub_cum_size(sub) + count(this%has_sub(1_i8:cell, sub), kind=i8)
  end function river_upscaler_cell_sub

  subroutine river_upscaler_destroy(this)
    implicit none
    class(river_upscaler_t), intent(inout) :: this

    if (allocated(this%scc_gauges)) deallocate(this%scc_gauges)
    if (allocated(this%scc_coarse_gauges)) deallocate(this%scc_coarse_gauges)
    if (allocated(this%scc_map)) deallocate(this%scc_map)
    if (allocated(this%is_scc_gauge)) deallocate(this%is_scc_gauge)
    if (allocated(this%is_scc_coarse_gauge)) deallocate(this%is_scc_coarse_gauge)
    if (allocated(this%leaving_cells)) deallocate(this%leaving_cells)
    if (allocated(this%has_sub)) deallocate(this%has_sub)
    if (allocated(this%sub_size)) deallocate(this%sub_size)
    if (allocated(this%n_sub_nodes)) deallocate(this%n_sub_nodes)
    if (allocated(this%node_sub)) deallocate(this%node_sub)
    if (allocated(this%sink_map)) deallocate(this%sink_map)
    if (allocated(this%link_start)) deallocate(this%link_start)
    if (allocated(this%is_link_start)) deallocate(this%is_link_start)
    if (allocated(this%stream_mask)) deallocate(this%stream_mask)
    if (allocated(this%stream_sub)) deallocate(this%stream_sub)
    if (allocated(this%link_end)) deallocate(this%link_end)

  end subroutine river_upscaler_destroy

end module mo_river_upscaler
