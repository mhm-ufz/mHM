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
  use mo_river, only: river_t
  use mo_grid, only: grid_t, cartesian
  use mo_grid_scaler, only: scaler_t, down_scaling
  use mo_message, only: error_message, message
  use mo_utils, only: optval

  implicit none
  private

  integer(i4), public, parameter :: upscale_legacy = 0_i4 !< connect to the first entered coarse cell
  integer(i4), public, parameter :: upscale_flow = 1_i4 !< connect to the next downstream link endpoint

  !> Temporary construction data, released when river_upscaler_t%init returns.
  type :: river_upscaler_scratch_t
    type(scaler_t) :: scaler
    integer(i4) :: nsub
    integer(i4), allocatable :: scc_map(:)
    integer(i8), allocatable :: scc_gauge_cells(:)
    logical, allocatable :: is_scc_gauge(:)
  end type river_upscaler_scratch_t

  !> \class river_upscaler_t
  !> \brief River network upscaler
  !> \details upscale river network respecting scc.
  !! Coarse nodes are ordered first by sub-catchment ID and then by increasing coarse-cell ID.
  type, public :: river_upscaler_t
    type(river_t), pointer :: fine_river => null() !< river definition at fine grid
    type(river_t), pointer :: coarse_river => null() !< river definition at coarse grid
    logical, allocatable :: stream_mask(:) !< mask marking the upscaled stream at fine grid size(fine\%ncells)
    integer(i8), allocatable :: link_start(:) !< starting cell at fine grid of coarse river node (0 for sink) size(coarse\%n_nodes)
    integer(i8), allocatable :: link_end(:) !< downstream fine endpoint of coarse river node (0 for sink) size(coarse\%n_nodes)
    integer(i8), allocatable :: scc_coarse_gauges(:) !< coarse node ID for each SCC gauge size(nsub-1)
  contains
    procedure, public :: init => river_upscaler_init
    procedure, private :: init_scc_gauges => river_upscaler_init_scc_gauges
    procedure, private :: upscale => river_upscaler_upscale
    procedure, private :: write_diagnostics => river_upscaler_write_diagnostics
    procedure, public :: calc_celerity => river_upscaler_celerity
    procedure, public :: destroy => river_upscaler_destroy
  end type river_upscaler_t

contains

  !> \brief Setup river upscaler from fine river and coarse target grid.
  subroutine river_upscaler_init(this, &
    fine_river, coarse_river, coarse_grid, scc_gauges, scc_latlon, upscale_mode, length_percentile, tol, &
    diagnostics_path, retain_stream_mask)
    implicit none
    class(river_upscaler_t), target, intent(inout) :: this
    type(river_t), pointer, intent(in) :: fine_river !< river definition at fine grid
    type(river_t), pointer, intent(in) :: coarse_river !< pointer to coarse river definition to be determined
    type(grid_t), pointer, intent(in) :: coarse_grid !< coarse grid for the upscaled river network
    real(dp), dimension(:,:), intent(in), optional :: scc_gauges !< gauge locations on fine river dim 1: id, dim 2: (x,y)
    logical, optional, intent(in) :: scc_latlon !< Whether the scc gauges coordinates are geographical (default: .true.)
    integer(i4), optional, intent(in) :: upscale_mode !< upscaling mode (0: legacy, 1: FLOW-like; default: 1)
    real(dp), optional, intent(in) :: length_percentile !< [%] percentile for lower cut-off of upscaled link-length (40 by default)
    real(dp), optional, intent(in) :: tol !< tolerance for cell factor comparison (default: 1.e-7)
    character(*), optional, intent(in) :: diagnostics_path !< optional diagnostics output file
    logical, optional, intent(in) :: retain_stream_mask !< retain fine stream mask for variable celerity
    type(river_upscaler_scratch_t) :: scratch
    integer(i4) :: mode
    real(dp) :: len_percentile
    logical :: keep_stream_mask

    mode = optval(upscale_mode, upscale_flow)
    len_percentile = optval(length_percentile, 40.0_dp)
    if (mode /= upscale_legacy .and. mode /= upscale_flow) then
      call error_message("river_upscaler: upscale_mode needs to be 0 (legacy) or 1 (FLOW-like)")
    end if
    if (len_percentile < 0.0_dp .or. len_percentile > 100.0_dp) then
      call error_message("river_upscaler: length_percentile needs to be in the range [0, 100]")
    end if

    if (allocated(this%link_start)) deallocate(this%link_start)
    if (allocated(this%link_end)) deallocate(this%link_end)
    if (allocated(this%stream_mask)) deallocate(this%stream_mask)
    if (allocated(this%scc_coarse_gauges)) deallocate(this%scc_coarse_gauges)
    this%fine_river => fine_river
    this%coarse_river => coarse_river
    this%coarse_river%grid => coarse_grid
    call scratch%scaler%init(this%fine_river%grid, this%coarse_river%grid, tol=tol)
    if (scratch%scaler%scaling_mode == down_scaling) then
      call error_message("river_upscaler: target grid needs to be coarser then input")
    end if
    ! TODO: shortcut for same resolution of fine and coarse river

    ! sanity check to have facc available for upscaling
    if (.not.allocated(this%fine_river%facc)) then
      if (.not.allocated(this%fine_river%order%id)) call this%fine_river%calc_order(root=.true.)
      call this%fine_river%calc_facc()
    end if

    ! initialize scc related variables
    call message("river_upscaler: initialize scc")
    call this%init_scc_gauges(scc_gauges, scc_latlon, scratch)

    ! upscale graph and calculate stream features in one fine-river trace
    call message("river_upscaler: upscale river graph and stream features")
    keep_stream_mask = optval(retain_stream_mask, .false.)
    call this%upscale(mode, len_percentile, scratch, diagnostics_path, keep_stream_mask)

  end subroutine river_upscaler_init

  !> \brief Initialize SCC related variables
  subroutine river_upscaler_init_scc_gauges(this, scc_gauges, scc_latlon, scratch)
    class(river_upscaler_t), target, intent(inout) :: this
    real(dp), dimension(:,:), intent(in), optional :: scc_gauges !< gauge locations on fine river dim 1: id, dim 2: (x,y)
    logical, optional, intent(in) :: scc_latlon !< Whether the scc gauges coordinates are geographical (default: .true.)
    type(river_upscaler_scratch_t), intent(inout) :: scratch
    integer(i4) :: i, n
    integer(i8) :: k
    logical :: aux, latlon

    if (.not.present(scc_gauges)) then
      call message("river_upscaler%init_scc: no scc gauges provided, initialize without scc")
      scratch%nsub = 1_i4
      allocate(scratch%scc_gauge_cells(0))
      allocate(scratch%scc_map(this%fine_river%n_nodes))
      allocate(scratch%is_scc_gauge(this%fine_river%n_nodes))
      ! set default
      !$omp parallel do default(shared)
      do k = 1_i8, this%fine_river%n_nodes
        scratch%scc_map(k) = 1_i4
        scratch%is_scc_gauge(k) = .false.
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
    scratch%nsub = n + 1_i4
    if (this%fine_river%scc) call error_message("river_upscaler%init_scc: need a D8 river to initialize SCC")
    if (.not.allocated(this%fine_river%facc)) call error_message("river_upscaler%init_scc: facc not available")

    ! find gauge cell ids
    allocate(scratch%scc_gauge_cells(n))
    call message("river_upscaler%init_scc: find gauge cell ids")
    if (aux) then
      call message("river_upscaler%init_scc: find gauge cell with aux coordinates (slow)")
      !$omp parallel do default(shared)
      do i = 1_i4, n
        scratch%scc_gauge_cells(i) = this%fine_river%grid%closest_cell_id(scc_gauges(i,:), use_aux=aux)
      end do
      !$omp end parallel do
    else
      scratch%scc_gauge_cells = this%fine_river%grid%closest_cell_id_by_axes(scc_gauges)
    end if

    ! calculate scc map
    if (.not.allocated(this%fine_river%order%id)) call this%fine_river%calc_order(root=.true.)
    call this%fine_river%label_subcatchments( &
      scratch%scc_map, scratch%scc_gauge_cells, default_label=scratch%nsub) ! base catchment gets id n+1

    ! determine gauges
    allocate(scratch%is_scc_gauge(this%fine_river%n_nodes))
    call message("river_upscaler%init_scc: determine scc gauges on fine grid")
    ! set default
    !$omp parallel do default(shared)
    do k = 1_i8, this%fine_river%n_nodes
      scratch%is_scc_gauge(k) = .false.
    end do
    !$omp end parallel do
    ! fill in gauges
    !$omp parallel do default(shared)
    do i = 1_i4, n
      scratch%is_scc_gauge(scratch%scc_gauge_cells(i)) = .true.
    end do
    !$omp end parallel do
    this%coarse_river%scc = scratch%nsub > 1_i4
  end subroutine river_upscaler_init_scc_gauges

  !> \brief Setup the coarse graph and stream features with optional SCC nodes.
  !> \details Both modes trace source-inclusive and endpoint-exclusive fine links.
  subroutine river_upscaler_upscale(this, upscale_mode, length_percentile, scratch, diagnostics_path, retain_stream_mask)
    use mo_percentile, only: percentile
    use mo_utils, only: prefix_sum
    implicit none
    class(river_upscaler_t), target, intent(inout) :: this
    integer(i4), intent(in) :: upscale_mode
    real(dp), intent(in) :: length_percentile
    type(river_upscaler_scratch_t), intent(inout) :: scratch
    character(*), optional, intent(in) :: diagnostics_path
    logical, intent(in) :: retain_stream_mask
    integer(i8), allocatable :: down(:), ids(:,:)
    integer(i4), allocatable :: facc(:,:), scc_map(:,:)
    integer(i8), allocatable :: node_from_cell_sub(:,:), sink_map(:)
    integer(i8), allocatable :: sub_size(:), sub_cum_size(:)
    integer(i4), allocatable :: node_sub(:)
    logical, allocatable :: leave_mask(:,:), has_sub(:,:), leaving_cells(:), is_link_start(:), stream_mask(:), &
      is_scc_coarse_gauge(:)
    integer(i8) :: i, k, node, next, cell, facc_max_i, n_nodes, n_links
    integer(i4) :: j, sub, facc_max, ix, iy, loc(2)
    integer(i4) :: yl, yu, xl, xu
    real(dp) :: length_cutoff
    real(dp), allocatable :: node_x(:), node_y(:)
    logical :: has_sub_init, mark_stream

    if (.not.allocated(this%fine_river%facc)) call error_message("river_upscaler%upscale: facc not available")
    if (.not.allocated(this%fine_river%link_length)) call error_message("river_upscaler%upscale: link length not available")
    if (this%fine_river%points%n_points /= this%fine_river%n_nodes) then
      call error_message("river_upscaler%upscale: node location not available")
    end if

    call message("river_upscaler%upscale: initialize sub-catchment attributes")
    allocate(has_sub(this%coarse_river%grid%ncells, scratch%nsub))
    has_sub_init = .not.this%coarse_river%scc
    !$omp parallel do default(shared)
    do j = 1_i4, scratch%nsub
      has_sub(:,j) = has_sub_init
    end do
    !$omp end parallel do

    call message("river_upscaler%upscale: find leaving cells")
    allocate(leaving_cells(this%fine_river%grid%ncells))
    !$omp parallel do default(shared)
    do i = 1_i8, this%fine_river%grid%ncells
      if (this%fine_river%is_sink(i)) then
        leaving_cells(i) = .false.
      else
        leaving_cells(i) = scratch%scaler%id_map(i) /= scratch%scaler%id_map(this%fine_river%down(i))
      end if
      if (this%coarse_river%scc) then
        !$omp atomic write
        has_sub(scratch%scaler%id_map(i), scratch%scc_map(i)) = .true.
      end if
    end do
    !$omp end parallel do

    call message("river_upscaler%upscale: allocate coarse river attributes")
    allocate(sub_size(scratch%nsub))
    allocate(sub_cum_size(scratch%nsub))
    !$omp parallel do default(shared)
    do j = 1_i4, scratch%nsub
      sub_size(j) = count(has_sub(:,j), kind=i8)
    end do
    !$omp end parallel do
    call prefix_sum(sub_size, sub_cum_size, shift=1_i8, start=0_i8)
    n_nodes = sub_cum_size(scratch%nsub) + sub_size(scratch%nsub) ! total number of coarse nodes

    call message("river_upscaler%upscale: create cell/sub to node map")
    allocate(node_from_cell_sub(this%coarse_river%grid%ncells, scratch%nsub))
    !$omp parallel do default(shared) private(i,k)
    do j = 1_i4, scratch%nsub
      k = sub_cum_size(j)
      do i = 1_i8, this%coarse_river%grid%ncells
        if (has_sub(i,j)) then
          k = k + 1_i8
          node_from_cell_sub(i,j) = k
        else
          node_from_cell_sub(i,j) = 0_i8
        end if
      end do
    end do
    !$omp end parallel do
    deallocate(has_sub, sub_size, sub_cum_size)

    ! coarse river attributes
    if (this%coarse_river%scc) then
      allocate(this%coarse_river%node_cell(n_nodes))
      allocate(this%coarse_river%area_fraction(n_nodes))
    end if
    allocate(this%coarse_river%is_sink(n_nodes))
    allocate(this%coarse_river%link_length(n_nodes))
    allocate(node_x(n_nodes), node_y(n_nodes))
    ! upscaler attributes
    allocate(is_scc_coarse_gauge(n_nodes))
    allocate(node_sub(n_nodes))
    allocate(sink_map(n_nodes))
    allocate(this%link_start(n_nodes))
    allocate(this%link_end(n_nodes))

    ! initialize attributes
    !$omp parallel do default(shared)
    do k = 1_i8, n_nodes
      this%coarse_river%link_length(k) = 0.0_dp
      this%coarse_river%is_sink(k) = .false.
      is_scc_coarse_gauge(k) = .false.
      sink_map(k) = 0_i8
      this%link_start(k) = 0_i8
      this%link_end(k) = 0_i8
    end do
    !$omp end parallel do

    ! determine sub-catchment for each coarse node
    call message("river_upscaler%upscale: determine sub-catchment for each coarse node")
    !$omp parallel do default(shared) private(k, i)
    do j = 1_i4, scratch%nsub
      do i = 1_i8, this%coarse_river%grid%ncells
        k = node_from_cell_sub(i, j)
        if (k == 0_i8) cycle
        if (this%coarse_river%scc) this%coarse_river%node_cell(k) = i
        node_sub(k) = j
      end do
    end do
    !$omp end parallel do

    ! determine coarse gauge for each sub-catchment gauge
    call message("river_upscaler%upscale: determine coarse gauge for each sub-catchment gauge")
    allocate(this%scc_coarse_gauges(size(scratch%scc_gauge_cells)))
    !$omp parallel do default(shared) private(k)
    do j = 1_i4, size(scratch%scc_gauge_cells)
      k = node_from_cell_sub(scratch%scaler%id_map(scratch%scc_gauge_cells(j)), j)
      this%scc_coarse_gauges(j) = k
      is_scc_coarse_gauge(k) = .true.
    end do
    !$omp end parallel do

    ! unpack scc map for fine river
    allocate(scc_map(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%unpack_into(scratch%scc_map, scc_map)

    ! determine area fraction of each coarse node
    if (this%coarse_river%scc) then
      call message("river_upscaler%upscale: determine area fraction of each coarse node")
      !$omp parallel do default(shared)
      do k = 1_i8, n_nodes
        this%coarse_river%area_fraction(k) = scratch%scaler%cell_fraction( &
          class_map   = scc_map, &
          coarse_cell = this%coarse_river%node_cell(k), &
          class_id    = node_sub(k))
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
      k = scratch%scaler%id_map(i)
      sub = scratch%scc_map(i)
      ! sinks in base catchment are only sinks if they are the maximum facc in the coarse cell
      if (sub == scratch%nsub) then
        call scratch%scaler%coarse_bounds(k, xl, xu, yl, yu)
        if (this%fine_river%facc(i) < maxval(facc(xl:xu,yl:yu), mask=(scc_map(xl:xu,yl:yu)==sub))) cycle
      end if
      node = node_from_cell_sub(k, sub)
      this%coarse_river%is_sink(node) = .true.
      sink_map(node) = i
    end do
    !$omp end parallel do

    ! unpack leaving cells for fine river
    allocate(leave_mask(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%unpack_into(leaving_cells, leave_mask)
    ! create ID matrix for fine river
    allocate(ids(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%gen_id_matrix(ids)

    ! determine starting fine cell for each coarse node
    call message("river_upscaler%upscale: determine starting fine cell for each coarse node")
    !$omp parallel do default(shared) private(k, sub, xl, xu, yl, yu, loc, ix, iy)
    do i = 1_i8, n_nodes
      if (this%coarse_river%is_sink(i)) cycle
      sub = node_sub(i)
      if (is_scc_coarse_gauge(i)) then
        this%link_start(i) = scratch%scc_gauge_cells(sub)
      else
        if (this%coarse_river%scc) then
          k = this%coarse_river%node_cell(i)
          call scratch%scaler%coarse_bounds(k, xl, xu, yl, yu)
          loc = maxloc(facc(xl:xu,yl:yu), mask=(leave_mask(xl:xu,yl:yu).and.(scc_map(xl:xu,yl:yu)==sub)))
        else
          k = i
          call scratch%scaler%coarse_bounds(k, xl, xu, yl, yu)
          loc = maxloc(facc(xl:xu,yl:yu), mask=(leave_mask(xl:xu,yl:yu)))
        end if
        ix = xl + loc(1) - 1_i4
        iy = yl + loc(2) - 1_i4
        this%link_start(i) = ids(ix, iy)
      end if
    end do
    !$omp end parallel do
    deallocate(ids, facc, scc_map, leave_mask)
    deallocate(scratch%scc_gauge_cells, is_scc_coarse_gauge, node_sub)
    if (.not.present(diagnostics_path)) deallocate(leaving_cells)

    ! initialize link starts and, when required, the fine stream mask
    mark_stream = retain_stream_mask .or. present(diagnostics_path)
    allocate(is_link_start(this%fine_river%n_nodes))
    if (mark_stream) allocate(stream_mask(this%fine_river%n_nodes))
    !$omp parallel do default(shared)
    do i = 1_i8, this%fine_river%n_nodes
      is_link_start(i) = .false.
      if (mark_stream) stream_mask(i) = .false.
    end do
    !$omp end parallel do

    ! initialize coarse node locations and mark link start cells
    call message("river_upscaler%upscale: initialize coarse node locations and mark link start cells")
    !$omp parallel do default(shared)
    do i = 1_i8, n_nodes
      if (this%coarse_river%is_sink(i)) then
        node_x(i) = this%fine_river%points%x(sink_map(i))
        node_y(i) = this%fine_river%points%y(sink_map(i))
      else
        node_x(i) = this%fine_river%points%x(this%link_start(i))
        node_y(i) = this%fine_river%points%y(this%link_start(i))
        is_link_start(this%link_start(i)) = .true.
      end if
    end do
    !$omp end parallel do
    call this%coarse_river%points%init(node_x, node_y, coordsys=this%fine_river%points%coordsys)

    ! mark stream mask and determine coarse link lengths and downstream nodes
    call message("river_upscaler%upscale: mark stream mask and determine coarse link lengths and downstream nodes")
    allocate(down(n_nodes))
    ! Link traces vary strongly in length; dynamic assignment avoids a long static-schedule tail.
    !$omp parallel do default(shared) private(cell, next, sub) schedule(dynamic, 1)
    do i = 1_i8, n_nodes
      if (this%coarse_river%is_sink(i)) then
        if (mark_stream) then
          !$omp atomic write
          stream_mask(sink_map(i)) = .true.
        end if
        down(i) = 0_i8
        cycle
      end if
      cell = this%link_start(i)
      next = this%fine_river%down(cell)
      if (upscale_mode == upscale_legacy) then
        sub = scratch%scc_map(next)
        down(i) = node_from_cell_sub(scratch%scaler%id_map(next), sub)
      end if
      do
        if (mark_stream) then
          !$omp atomic write
          stream_mask(cell) = .true.
        end if
        this%coarse_river%link_length(i) = this%coarse_river%link_length(i) + this%fine_river%link_length(cell)
        cell = this%fine_river%down(cell)
        if (is_link_start(cell)) exit
        if (scratch%is_scc_gauge(cell)) exit
        if (this%fine_river%is_sink(cell)) exit
      end do
      if (mark_stream) then
        !$omp atomic write
        stream_mask(cell) = .true.
      end if
      this%link_end(i) = cell
      if (upscale_mode == upscale_flow) then
        sub = scratch%scc_map(cell)
        down(i) = node_from_cell_sub(scratch%scaler%id_map(cell), sub)
      end if
    end do
    !$omp end parallel do
    deallocate(scratch%is_scc_gauge)
    if (.not.present(diagnostics_path)) deallocate(scratch%scc_map, is_link_start)

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
        do j = 1_i4, scratch%nsub
          k = node_from_cell_sub(i, j)
          if (k == 0_i8) cycle
          if (this%coarse_river%is_sink(k)) then
            node = sink_map(k)
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
    deallocate(node_from_cell_sub, sink_map)

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

    if (present(diagnostics_path)) then
      call this%write_diagnostics( &
        diagnostics_path, scratch%scaler, scratch%scc_map, leaving_cells, stream_mask, is_link_start)
    end if
    if (retain_stream_mask) call move_alloc(stream_mask, this%stream_mask)

  end subroutine river_upscaler_upscale

  !> \brief Write river-upscaling diagnostics while construction data is available.
  subroutine river_upscaler_write_diagnostics(this, path, scaler, sub_map, leaving_cells, stream_mask, is_link_start)
    implicit none
    class(river_upscaler_t), target, intent(in) :: this
    character(*), intent(in) :: path
    type(scaler_t), intent(in) :: scaler
    integer(i4), intent(in) :: sub_map(:)
    logical, intent(in) :: leaving_cells(:), stream_mask(:), is_link_start(:)
    integer(i4), allocatable :: stream_sub(:)
    logical, allocatable :: highlight(:)
    integer(i8) :: i

    call message("river_upscaler: write diagnostics to file: " // trim(path))
    allocate(stream_sub(this%fine_river%n_nodes), source=nodata_i4)
    allocate(highlight(this%fine_river%n_nodes))
    !$omp parallel do default(shared)
    do i = 1_i8, this%fine_river%n_nodes
      if (stream_mask(i)) stream_sub(i) = sub_map(i)
      highlight(i) = is_link_start(i) .or. this%fine_river%is_sink(i)
    end do
    !$omp end parallel do
    call this%fine_river%export( &
      path        = path, &
      sub_map     = sub_map, &
      leaving     = leaving_cells, &
      stream_mask = stream_mask, &
      stream_sub  = stream_sub, &
      highlight   = highlight, &
      factor      = scaler%factor &
    )
  end subroutine river_upscaler_write_diagnostics

  !> \brief calculate the celerity c_i from slope s_i (i - cell index)
  subroutine river_upscaler_celerity(this, gamma, constant_celerity, slope)
    implicit none
    class(river_upscaler_t), target, intent(inout) :: this
    real(dp), intent(in) :: gamma !< model parameter: c_i = gamma * sqrt(s_i) or c = gamma
    logical, optional, intent(in) :: constant_celerity !< whether celerity is assumed constant: c = gamma (default: .false.)
    real(dp), optional, intent(in) :: slope(:) !< [%] river slope on fine grid: size(fine\%ncells)
    integer(i8) :: i, cell
    real(dp) :: n

    if (optval(constant_celerity, .false.)) then
      call message("river_upscaler: constant celerity assumed, set c_i = gamma for all coarse nodes")
      ! constant celerity, no need to calculate from fine river
      call this%coarse_river%calc_celerity(gamma, constant_celerity)
      return
    end if

    if (.not.allocated(this%stream_mask)) then
      call error_message("river_upscaler%calc_celerity: variable celerity requires retain_stream_mask=.true. during init")
    end if

    ! first calculate celerity on fine river
    call this%fine_river%calc_celerity(gamma, constant_celerity, slope, this%stream_mask)

    call message("river_upscaler: calculate celerity on coarse river from fine river")
    if (.not.allocated(this%coarse_river%celerity)) allocate(this%coarse_river%celerity(this%coarse_river%n_nodes))
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

  subroutine river_upscaler_destroy(this)
    implicit none
    class(river_upscaler_t), target, intent(inout) :: this

    if (allocated(this%link_start)) deallocate(this%link_start)
    if (allocated(this%stream_mask)) deallocate(this%stream_mask)
    if (allocated(this%link_end)) deallocate(this%link_end)
    if (allocated(this%scc_coarse_gauges)) deallocate(this%scc_coarse_gauges)

  end subroutine river_upscaler_destroy

end module mo_river_upscaler
