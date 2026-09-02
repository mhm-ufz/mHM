!> \file    mo_river_upscaler.f90
!> \copydoc mo_river_upscaler

!> \brief   River upscaler.
!> \details This module contains an upscaler for river networks.
!> \version 0.1
!> \changelog
!! - Luis Samaniego (2005/2012): original routing-network upscaling.
!! - Rohini Kumar (2014): grid-geometry support for upscaling.
!! - Stephan Thober (2015-2020): mRM port, multi-outlet support, and link-length cutoff.
!! - Robert Schweppe (2018): prior routing-network refactoring.
!! - Pallav Shrestha (2018-2023): lake-aware SCC topology and scalability extensions.
!! - Sebastian Mueller (2025-2026): v6 sparse DAG-based upscaler rewrite.
!> \authors Luis Samaniego, Rohini Kumar, Stephan Thober, Robert Schweppe, Sebastian Mueller, Pallav Shrestha
!> \date    May 2025
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_upscaler

  use mo_kind, only: i4, i8, dp
  use mo_constants, only: nodata_i4
  use mo_river, only: river_t
  use mo_grid, only: grid_t
  use mo_grid_scaler, only: scaler_t, down_scaling, no_scaling
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
    integer(i4) :: n_lakes = 0_i4
    integer(i4), allocatable :: scc_map(:)
    integer(i8), allocatable :: scc_gauge_cells(:)
    integer(i8), allocatable :: lake_ids(:)
    logical, allocatable :: is_scc_gauge(:)
    logical, allocatable :: is_lake_outlet(:)
  end type river_upscaler_scratch_t

  !> Sparse mapping from an existing (coarse cell, SCC label) pair to its coarse river node.
  !> \authors Sebastian Mueller, Pallav Shrestha
  type :: sparse_cell_sub_t
    integer(i8), allocatable :: offset(:) !< one-based CSR offsets size(ncells+1)
    integer(i4), allocatable :: sub(:) !< sorted SCC labels for all existing pairs
    integer(i8), allocatable :: node(:) !< coarse node for each pair
  contains
    procedure :: find => sparse_cell_sub_find
  end type sparse_cell_sub_t

  !> \class river_upscaler_t
  !> \brief River network upscaler
  !> \authors Sebastian Mueller, Pallav Shrestha
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
  !> \authors Sebastian Mueller, Pallav Shrestha
  subroutine river_upscaler_init(this, &
    fine_river, coarse_river, coarse_grid, scc_nodes, lake_outlet_nodes, lake_ids, upscale_mode, length_percentile, tol, &
    diagnostics_path, retain_stream_mask)
    implicit none
    class(river_upscaler_t), target, intent(inout) :: this !< Upscaler instance to initialize.
    type(river_t), pointer, intent(in) :: fine_river !< river definition at fine grid
    type(river_t), pointer, intent(in) :: coarse_river !< pointer to coarse river definition to be determined
    type(grid_t), pointer, intent(in) :: coarse_grid !< coarse grid for the upscaled river network
    integer(i8), intent(in), optional :: scc_nodes(:) !< additional SCC nodes already snapped to the fine river
    integer(i8), intent(in), optional :: lake_outlet_nodes(:) !< lake SCC nodes, ordered before additional SCC nodes
    integer(i8), intent(in), optional :: lake_ids(:) !< stable IDs aligned with lake_outlet_nodes
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
    call this%init_scc_gauges(scc_nodes, lake_outlet_nodes, lake_ids, scratch)

    ! upscale graph and calculate stream features in one fine-river trace
    call message("river_upscaler: upscale river graph and stream features")
    keep_stream_mask = optval(retain_stream_mask, .false.)
    call this%upscale(mode, len_percentile, scratch, diagnostics_path, keep_stream_mask)

  end subroutine river_upscaler_init

  !> \brief Initialize SCC related variables
  !> \authors Sebastian Mueller, Pallav Shrestha
  subroutine river_upscaler_init_scc_gauges(this, scc_nodes, lake_outlet_nodes, lake_ids, scratch)
    class(river_upscaler_t), target, intent(inout) :: this !< Upscaler holding the fine and coarse rivers.
    integer(i8), intent(in), optional :: scc_nodes(:) !< Snapped ordinary L0 SCC outlet nodes.
    integer(i8), intent(in), optional :: lake_outlet_nodes(:) !< Snapped L0 lake outlets, ordered before ordinary SCC nodes.
    integer(i8), intent(in), optional :: lake_ids(:) !< Stable lake IDs aligned with lake_outlet_nodes.
    type(river_upscaler_scratch_t), intent(inout) :: scratch !< Temporary SCC construction state.
    integer(i4) :: i, n, n_scc
    integer(i8) :: k
    logical :: has_scc, has_lakes

    has_scc = present(scc_nodes)
    if (has_scc) has_scc = size(scc_nodes, kind=i8) > 0_i8
    has_lakes = present(lake_outlet_nodes) .or. present(lake_ids)
    if (present(lake_outlet_nodes) .neqv. present(lake_ids)) &
      call error_message("river_upscaler%init_scc: lake outlet nodes and IDs must be provided together")
    if (has_lakes) then
      if (size(lake_outlet_nodes, kind=i8) /= size(lake_ids, kind=i8)) &
        call error_message("river_upscaler%init_scc: lake outlet node and ID counts differ")
      has_lakes = size(lake_ids, kind=i8) > 0_i8
    end if

    if (.not.has_scc .and. .not.has_lakes) then
      call message("river_upscaler%init_scc: no scc gauges provided, initialize without scc")
      scratch%nsub = 1_i4
      allocate(scratch%scc_gauge_cells(0))
      allocate(scratch%lake_ids(0))
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

    scratch%n_lakes = 0_i4
    if (has_lakes) scratch%n_lakes = size(lake_ids, kind=i4)
    n_scc = 0_i4
    if (has_scc) n_scc = size(scc_nodes, kind=i4)
    n = scratch%n_lakes + n_scc
    scratch%nsub = n + 1_i4
    if (this%fine_river%scc) call error_message("river_upscaler%init_scc: need a D8 river to initialize SCC")
    if (.not.allocated(this%fine_river%facc)) call error_message("river_upscaler%init_scc: facc not available")

    allocate(scratch%scc_gauge_cells(n))
    if (has_lakes) then
      scratch%scc_gauge_cells(:scratch%n_lakes) = lake_outlet_nodes
      scratch%lake_ids = lake_ids
    else
      allocate(scratch%lake_ids(0))
    end if
    if (has_scc) scratch%scc_gauge_cells(scratch%n_lakes + 1_i4:) = scc_nodes
    if (any(scratch%scc_gauge_cells < 1_i8) .or. any(scratch%scc_gauge_cells > this%fine_river%n_nodes)) &
      call error_message("river_upscaler%init_scc: SCC node outside fine river")

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
    if (scratch%n_lakes > 0_i4) then
      allocate(scratch%is_lake_outlet(this%fine_river%n_nodes))
      !$omp parallel do default(shared) schedule(static)
      do k = 1_i8, this%fine_river%n_nodes
        scratch%is_lake_outlet(k) = .false.
      end do
      !$omp end parallel do
      !$omp parallel do default(shared) schedule(static)
      do i = 1_i4, scratch%n_lakes
        scratch%is_lake_outlet(scratch%scc_gauge_cells(i)) = .true.
      end do
      !$omp end parallel do
    end if
    this%coarse_river%scc = scratch%nsub > 1_i4
  end subroutine river_upscaler_init_scc_gauges

  !> \brief Return the coarse node for an existing sparse (cell, SCC label) pair, or zero.
  !> \authors Sebastian Mueller, Pallav Shrestha
  pure integer(i8) function sparse_cell_sub_find(this, cell, sub) result(node)
    class(sparse_cell_sub_t), intent(in) :: this !< Sparse cell/subcatchment mapping to query.
    integer(i8), intent(in) :: cell !< Coarse-grid cell ID.
    integer(i4), intent(in) :: sub !< SCC subcatchment label.
    integer(i8) :: left, right, middle

    node = 0_i8
    if (cell < 1_i8 .or. cell >= size(this%offset, kind=i8)) return
    left = this%offset(cell)
    right = this%offset(cell + 1_i8) - 1_i8
    do while (left <= right)
      middle = left + (right - left) / 2_i8
      if (this%sub(middle) < sub) then
        left = middle + 1_i8
      else if (this%sub(middle) > sub) then
        right = middle - 1_i8
      else
        node = this%node(middle)
        return
      end if
    end do
  end function sparse_cell_sub_find

  !> \brief Build sparse existing (coarse cell, SCC label) pairs without dense cell-by-catchment storage.
  !> \authors Sebastian Mueller, Pallav Shrestha
  subroutine build_sparse_cell_sub(scaler, pair_map, nsub, pairs)
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    use mo_orderpack, only: sort
    use mo_utils, only: prefix_sum
    type(scaler_t), intent(in) :: scaler !< Fine-to-coarse grid scaler.
    integer(i4), intent(in) :: pair_map(:,:) !< Fine-grid SCC labels; negative values are excluded lake cells.
    integer(i4), intent(in) :: nsub !< Number of SCC labels including the base catchment.
    type(sparse_cell_sub_t), intent(out) :: pairs !< Deterministic sparse coarse-cell/SCC mapping.
    integer(i8), allocatable :: count_cell(:), seen(:,:), sub_count(:), sub_start(:), sub_next(:)
    integer(i8) :: cell, k, pos, n_pairs
    integer(i4) :: x, y, xl, xu, yl, yu, sub, thread, n_threads

    n_threads = 1_i4
    !$ n_threads = omp_get_max_threads()
    allocate(count_cell(scaler%coarse_grid%ncells))
    allocate(seen(nsub, n_threads))
    !$omp parallel do default(shared) schedule(static)
    do thread = 1_i4, n_threads
      seen(:, thread) = 0_i8
    end do
    !$omp end parallel do

    !$omp parallel default(shared) private(thread,cell,x,y,xl,xu,yl,yu,sub)
    thread = 1_i4
    !$ thread = omp_get_thread_num() + 1_i4
    !$omp do schedule(static)
    do cell = 1_i8, scaler%coarse_grid%ncells
      count_cell(cell) = 0_i8
      call scaler%coarse_bounds(cell, xl, xu, yl, yu)
      do y = yl, yu
        do x = xl, xu
          sub = pair_map(x, y)
          if (sub < 1_i4) cycle
          if (seen(sub, thread) == cell) cycle
          seen(sub, thread) = cell
          count_cell(cell) = count_cell(cell) + 1_i8
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

    allocate(pairs%offset(scaler%coarse_grid%ncells + 1_i8))
    call prefix_sum(count_cell, pairs%offset(:scaler%coarse_grid%ncells), shift=1_i8, start=0_i8)
    pairs%offset(:scaler%coarse_grid%ncells) = pairs%offset(:scaler%coarse_grid%ncells) + 1_i8
    pairs%offset(scaler%coarse_grid%ncells + 1_i8) = &
      pairs%offset(scaler%coarse_grid%ncells) + count_cell(scaler%coarse_grid%ncells)
    n_pairs = pairs%offset(scaler%coarse_grid%ncells + 1_i8) - 1_i8
    allocate(pairs%sub(n_pairs), pairs%node(n_pairs))

    !$omp parallel default(shared) private(thread,cell,pos,x,y,xl,xu,yl,yu,sub)
    thread = 1_i4
    !$ thread = omp_get_thread_num() + 1_i4
    !$omp do schedule(static)
    do cell = 1_i8, scaler%coarse_grid%ncells
      pos = pairs%offset(cell)
      call scaler%coarse_bounds(cell, xl, xu, yl, yu)
      do y = yl, yu
        do x = xl, xu
          sub = pair_map(x, y)
          if (sub < 1_i4) cycle
          if (seen(sub, thread) == -cell) cycle
          seen(sub, thread) = -cell
          pairs%sub(pos) = sub
          pos = pos + 1_i8
        end do
      end do
      if (pairs%offset(cell + 1_i8) - pairs%offset(cell) > 1_i8) &
        call sort(pairs%sub(pairs%offset(cell):pairs%offset(cell + 1_i8) - 1_i8))
    end do
    !$omp end do
    !$omp end parallel
    deallocate(seen, count_cell)

    allocate(sub_count(nsub))
    !$omp parallel do default(shared) schedule(static)
    do sub = 1_i4, nsub
      sub_count(sub) = 0_i8
    end do
    !$omp end parallel do
    !$omp parallel do default(shared) private(sub) schedule(static)
    do k = 1_i8, n_pairs
      sub = pairs%sub(k)
      !$omp atomic update
      sub_count(sub) = sub_count(sub) + 1_i8
    end do
    !$omp end parallel do
    allocate(sub_start(nsub), sub_next(nsub))
    call prefix_sum(sub_count, sub_start, shift=1_i8, start=0_i8)
    sub_next = sub_start
    ! Cell-major CSR traversal assigns increasing cell IDs within each subcatchment.
    do cell = 1_i8, scaler%coarse_grid%ncells
      do k = pairs%offset(cell), pairs%offset(cell + 1_i8) - 1_i8
        sub = pairs%sub(k)
        sub_next(sub) = sub_next(sub) + 1_i8
        pairs%node(k) = sub_next(sub)
      end do
    end do
    deallocate(sub_count, sub_start, sub_next)
  end subroutine build_sparse_cell_sub

  !> \brief Setup the coarse graph and stream features with optional SCC nodes.
  !> \authors Sebastian Mueller, Pallav Shrestha
  !> \details Both modes trace source-inclusive and endpoint-exclusive fine links.
  subroutine river_upscaler_upscale(this, upscale_mode, length_percentile, scratch, diagnostics_path, retain_stream_mask)
    use mo_percentile, only: percentile
    implicit none
    class(river_upscaler_t), target, intent(inout) :: this !< Initialized upscaler and rivers to populate.
    integer(i4), intent(in) :: upscale_mode !< Link-targeting mode: legacy or FLOW-like.
    real(dp), intent(in) :: length_percentile !< Lower percentile used to clip coarse link lengths.
    type(river_upscaler_scratch_t), intent(inout) :: scratch !< Temporary SCC and scaling construction state.
    character(*), optional, intent(in) :: diagnostics_path !< Optional fine-grid topology diagnostics path.
    logical, intent(in) :: retain_stream_mask !< Whether to retain the fine stream mask for celerity calculation.
    integer(i8), allocatable :: down(:), ids(:,:)
    integer(i4), allocatable :: facc(:,:), scc_map(:,:)
    integer(i8), allocatable :: sink_map(:)
    integer(i4), allocatable :: node_sub(:)
    logical, allocatable :: leave_mask(:,:), leaving_cells(:), is_link_start(:), stream_mask(:), &
      is_scc_coarse_gauge(:)
    type(sparse_cell_sub_t) :: sub_list !< Existing coarse-cell/SCC-label entries and their L3 nodes.
    integer(i8) :: i, k, p, node, next, cell, facc_max_i, n_nodes, n_links
    integer(i4) :: j, sub, lake_index, facc_max, ix, iy, loc(2)
    integer(i4) :: yl, yu, xl, xu
    real(dp) :: length_cutoff
    real(dp), allocatable :: node_x(:), node_y(:)
    logical :: mark_stream, lake_only_candidates, invalid_lake_label

    if (.not.allocated(this%fine_river%facc)) call error_message("river_upscaler%upscale: facc not available")
    if (.not.allocated(this%fine_river%link_length)) call error_message("river_upscaler%upscale: link length not available")
    if (this%fine_river%points%n_points /= this%fine_river%n_nodes) then
      call error_message("river_upscaler%upscale: node location not available")
    end if

    ! Pair-map convention during construction: positive SCC labels are eligible routing pairs;
    ! negative labels identify lake cells excluded from normal node construction.
    allocate(scc_map(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%unpack_into(scratch%scc_map, scc_map)
    if (scratch%n_lakes > 0_i4) then
      if (.not.allocated(this%fine_river%lake_map)) &
        call error_message("river_upscaler%upscale: lake outlets require a fine-river lake map")
      invalid_lake_label = .false.
      !$omp parallel do default(shared) private(ix,iy,sub) reduction(.or.:invalid_lake_label) schedule(static)
      do i = 1_i8, this%fine_river%n_nodes
        if (this%fine_river%lake_map(i) <= 0_i8) cycle
        sub = scratch%scc_map(i)
        if (sub < 1_i4 .or. sub > scratch%n_lakes) then
          invalid_lake_label = .true.
          cycle
        end if
        if (scratch%lake_ids(sub) /= this%fine_river%lake_map(i)) then
          invalid_lake_label = .true.
          cycle
        end if
        if (scratch%is_lake_outlet(i)) cycle
        ix = this%fine_river%grid%cell_ij(i, 1)
        iy = this%fine_river%grid%cell_ij(i, 2)
        scc_map(ix, iy) = -sub
      end do
      !$omp end parallel do
      if (invalid_lake_label) &
        call error_message("river_upscaler%upscale: lake map is inconsistent with lake SCC catchments")
    end if

    if (this%coarse_river%scc) then
      call message("river_upscaler%upscale: build sparse cell/sub-catchment pairs")
      call build_sparse_cell_sub(scratch%scaler, scc_map, scratch%nsub, sub_list)
      n_nodes = size(sub_list%sub, kind=i8)
    else
      ! Without SCCs every coarse cell has exactly one routing node. Construct
      ! this mapping directly instead of scanning the complete fine grid twice.
      call message("river_upscaler%upscale: initialize one node per coarse cell")
      n_nodes = this%coarse_river%grid%ncells
      allocate(sub_list%offset(n_nodes + 1_i8))
      allocate(sub_list%sub(n_nodes))
      allocate(sub_list%node(n_nodes))
      !$omp parallel do default(shared) schedule(static)
      do i = 1_i8, n_nodes
        sub_list%offset(i) = i
        sub_list%sub(i) = 1_i4
        sub_list%node(i) = i
      end do
      !$omp end parallel do
      sub_list%offset(n_nodes + 1_i8) = n_nodes + 1_i8
    end if

    call message("river_upscaler%upscale: find leaving cells")
    allocate(leaving_cells(this%fine_river%grid%ncells))
    !$omp parallel do default(shared)
    do i = 1_i8, this%fine_river%grid%ncells
      if (this%fine_river%is_sink(i)) then
        leaving_cells(i) = .false.
      else if (scratch%n_lakes > 0_i4) then
        if (this%fine_river%lake_map(i) > 0_i8 .and. .not.scratch%is_lake_outlet(i)) then
          leaving_cells(i) = .false.
        else
          leaving_cells(i) = scratch%scaler%id_map(i) /= scratch%scaler%id_map(this%fine_river%down(i)) .or. &
            this%fine_river%lake_map(this%fine_river%down(i)) > 0_i8
        end if
      else
        leaving_cells(i) = scratch%scaler%id_map(i) /= scratch%scaler%id_map(this%fine_river%down(i))
      end if
    end do
    !$omp end parallel do

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
    !$omp parallel do default(shared) private(k,p) schedule(static)
    do i = 1_i8, this%coarse_river%grid%ncells
      do p = sub_list%offset(i), sub_list%offset(i + 1_i8) - 1_i8
        k = sub_list%node(p)
        if (this%coarse_river%scc) this%coarse_river%node_cell(k) = i
        node_sub(k) = sub_list%sub(p)
      end do
    end do
    !$omp end parallel do

    ! determine coarse gauge for each sub-catchment gauge
    call message("river_upscaler%upscale: determine coarse gauge for each sub-catchment gauge")
    allocate(this%scc_coarse_gauges(size(scratch%scc_gauge_cells)))
    !$omp parallel do default(shared) private(k)
    do j = 1_i4, size(scratch%scc_gauge_cells)
      k = sub_list%find(scratch%scaler%id_map(scratch%scc_gauge_cells(j)), j)
      this%scc_coarse_gauges(j) = k
      if (k > 0_i8) is_scc_coarse_gauge(k) = .true.
    end do
    !$omp end parallel do
    if (any(this%scc_coarse_gauges == 0_i8)) &
      call error_message("river_upscaler%upscale: SCC outlet has no coarse node")

    ! Exclude every lake-surface cell from local runoff fractions and normal link-start selection.
    if (scratch%n_lakes > 0_i4) then
      !$omp parallel do default(shared) private(ix,iy) schedule(static)
      do i = 1_i8, this%fine_river%n_nodes
        if (this%fine_river%lake_map(i) <= 0_i8) cycle
        ix = this%fine_river%grid%cell_ij(i, 1)
        iy = this%fine_river%grid%cell_ij(i, 2)
        scc_map(ix, iy) = -scratch%scc_map(i)
      end do
      !$omp end parallel do
    end if

    ! determine area fraction of each coarse node
    if (this%coarse_river%scc) then
      call message("river_upscaler%upscale: determine area fraction of each coarse node")
      if (scratch%scaler%scaling_mode == no_scaling) then
        !$omp parallel do default(shared) private(xl,xu,yl,yu)
        do k = 1_i8, n_nodes
          call scratch%scaler%coarse_bounds(this%coarse_river%node_cell(k), xl, xu, yl, yu)
          this%coarse_river%area_fraction(k) = real(count( &
            (scc_map(xl:xu,yl:yu) == node_sub(k)) .and. this%fine_river%grid%mask(xl:xu,yl:yu)), dp) / &
            real(count(this%fine_river%grid%mask(xl:xu,yl:yu)), dp)
        end do
        !$omp end parallel do
      else
        !$omp parallel do default(shared)
        do k = 1_i8, n_nodes
          this%coarse_river%area_fraction(k) = scratch%scaler%cell_fraction( &
            class_map   = scc_map, &
            coarse_cell = this%coarse_river%node_cell(k), &
            class_id    = node_sub(k))
        end do
        !$omp end parallel do
      end if
      if (scratch%n_lakes > 0_i4) then
        allocate(this%coarse_river%cell_land_fraction(this%coarse_river%grid%ncells))
        !$omp parallel do default(shared) private(p, node) schedule(static)
        do cell = 1_i8, this%coarse_river%grid%ncells
          this%coarse_river%cell_land_fraction(cell) = 0.0_dp
          do p = sub_list%offset(cell), sub_list%offset(cell + 1_i8) - 1_i8
            node = sub_list%node(p)
            this%coarse_river%cell_land_fraction(cell) = this%coarse_river%cell_land_fraction(cell) + &
              this%coarse_river%area_fraction(node)
          end do
          if (this%coarse_river%cell_land_fraction(cell) > 0.0_dp) then
            do p = sub_list%offset(cell), sub_list%offset(cell + 1_i8) - 1_i8
              node = sub_list%node(p)
              this%coarse_river%area_fraction(node) = this%coarse_river%area_fraction(node) / &
                this%coarse_river%cell_land_fraction(cell)
            end do
          end if
        end do
        !$omp end parallel do
      end if
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
      node = sub_list%find(k, sub)
      if (node == 0_i8) cycle
      this%coarse_river%is_sink(node) = .true.
      sink_map(node) = i
    end do
    !$omp end parallel do

    ! unpack leaving cells for fine river
    allocate(leave_mask(this%fine_river%grid%nx, this%fine_river%grid%ny))
    call this%fine_river%grid%unpack_into(leaving_cells, leave_mask)
    ! Avoid a masked prefix-count lookup for every selected coarse-node start.
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
    deallocate(facc, leave_mask, ids)
    deallocate(is_scc_coarse_gauge, node_sub)
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
    this%coarse_river%points%coordsys = this%fine_river%points%coordsys
    this%coarse_river%points%n_points = n_nodes
    call move_alloc(node_x, this%coarse_river%points%x)
    call move_alloc(node_y, this%coarse_river%points%y)

    ! mark stream mask and determine coarse link lengths and downstream nodes
    call message("river_upscaler%upscale: mark stream mask and determine coarse link lengths and downstream nodes")
    allocate(down(n_nodes))
    ! Link traces vary strongly in length; dynamic assignment avoids a long static-schedule tail.
    !$omp parallel do default(shared) private(cell,next,sub) schedule(dynamic, 1)
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
        if (scratch%n_lakes > 0_i4) then
          if (this%fine_river%lake_map(next) > 0_i8) then
            down(i) = this%scc_coarse_gauges(sub)
          else
            down(i) = sub_list%find(scratch%scaler%id_map(next), sub)
          end if
        else
          down(i) = sub_list%find(scratch%scaler%id_map(next), sub)
        end if
      end if
      do
        if (mark_stream) then
          !$omp atomic write
          stream_mask(cell) = .true.
        end if
        this%coarse_river%link_length(i) = this%coarse_river%link_length(i) + this%fine_river%link_length(cell)
        cell = this%fine_river%down(cell)
        if (scratch%n_lakes > 0_i4) then
          if (this%fine_river%lake_map(cell) > 0_i8) exit
        end if
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
        if (scratch%n_lakes > 0_i4) then
          if (this%fine_river%lake_map(cell) > 0_i8) then
            down(i) = this%scc_coarse_gauges(sub)
          else
            down(i) = sub_list%find(scratch%scaler%id_map(cell), sub)
          end if
        else
          down(i) = sub_list%find(scratch%scaler%id_map(cell), sub)
        end if
      end if
    end do
    !$omp end parallel do
    if (any(down == 0_i8 .and. .not.this%coarse_river%is_sink)) then
      call error_message("river_upscaler: a non-sink routing node has no downstream node")
    end if
    deallocate(scratch%is_scc_gauge)
    if (.not.present(diagnostics_path)) deallocate(scratch%scc_map, is_link_start)

    if (allocated(this%coarse_river%fdir)) deallocate(this%coarse_river%fdir)
    call message("river_upscaler: initialize coarse river branching DAG")
    call this%coarse_river%init(down)
    deallocate(down)
    if (scratch%n_lakes > 0_i4) &
      call this%coarse_river%set_lake_nodes(this%scc_coarse_gauges(:scratch%n_lakes), scratch%lake_ids)

    if (this%coarse_river%scc) then
      call message("river_upscaler: find representative node for each coarse cell")
      allocate(this%coarse_river%cell_node_select(this%coarse_river%grid%ncells))
      !$omp parallel do default(shared) private(p,k,node,lake_index,facc_max,facc_max_i,xl,xu,yl,yu,ix,iy,lake_only_candidates)
      do i = 1_i8, this%coarse_river%grid%ncells
        facc_max = 0_i4
        facc_max_i = 0_i8
        lake_only_candidates = scratch%n_lakes > 0_i4
        if (lake_only_candidates) then
          do p = sub_list%offset(i), sub_list%offset(i + 1_i8) - 1_i8
            if (this%coarse_river%lake_id(sub_list%node(p)) == 0_i8) then
              lake_only_candidates = .false.
              exit
            end if
          end do
        end if
        do p = sub_list%offset(i), sub_list%offset(i + 1_i8) - 1_i8
          k = sub_list%node(p)
          if (this%coarse_river%is_sink(k)) then
            node = sink_map(k)
          else
            node = this%link_start(k)
          end if
          if (this%fine_river%facc(node) > facc_max) then
            facc_max = this%fine_river%facc(node)
            facc_max_i = k
          else if (lake_only_candidates .and. this%fine_river%facc(node) == facc_max) then
            if (facc_max_i == 0_i8) then
              facc_max_i = k
            else if (this%coarse_river%lake_id(k) < this%coarse_river%lake_id(facc_max_i)) then
              facc_max_i = k
            end if
          end if
        end do
        if (facc_max_i == 0_i8 .and. scratch%n_lakes > 0_i4) then
          call scratch%scaler%coarse_bounds(i, xl, xu, yl, yu)
          do iy = yl, yu
            do ix = xl, xu
              if (scc_map(ix, iy) > -1_i4 .or. scc_map(ix, iy) < -scratch%n_lakes) cycle
              lake_index = -scc_map(ix, iy)
              node = scratch%scc_gauge_cells(lake_index)
              if (this%fine_river%facc(node) > facc_max) then
                facc_max = this%fine_river%facc(node)
                facc_max_i = this%scc_coarse_gauges(lake_index)
              else if (this%fine_river%facc(node) == facc_max) then
                if (facc_max_i == 0_i8) then
                  facc_max_i = this%scc_coarse_gauges(lake_index)
                else if (scratch%lake_ids(lake_index) < this%coarse_river%lake_id(facc_max_i)) then
                  facc_max_i = this%scc_coarse_gauges(lake_index)
                end if
              end if
            end do
          end do
        end if
        this%coarse_river%cell_node_select(i) = facc_max_i
      end do
      !$omp end parallel do
      if (any(this%coarse_river%cell_node_select == 0_i8)) &
        call error_message("river_upscaler: coarse cell has no representative river node")
    end if
    deallocate(sink_map, scc_map)
    if (allocated(sub_list%offset)) deallocate(sub_list%offset, sub_list%sub, sub_list%node)

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
        diagnostics_path, scratch%scaler, scratch%scc_map, leaving_cells, stream_mask, is_link_start, &
        scratch%scc_gauge_cells(:scratch%n_lakes), scratch%lake_ids)
    end if
    if (retain_stream_mask) call move_alloc(stream_mask, this%stream_mask)
    deallocate(scratch%scc_gauge_cells)
    if (allocated(scratch%is_lake_outlet)) deallocate(scratch%is_lake_outlet)
    if (allocated(scratch%lake_ids)) deallocate(scratch%lake_ids)

  end subroutine river_upscaler_upscale

  !> \brief Write river-upscaling diagnostics while construction data is available.
  !> \authors Sebastian Mueller, Pallav Shrestha
  subroutine river_upscaler_write_diagnostics( &
    this, path, scaler, sub_map, leaving_cells, stream_mask, is_link_start, lake_outlet_cells, lake_ids)
    implicit none
    class(river_upscaler_t), target, intent(in) :: this !< Upscaler containing the fine river to export.
    character(*), intent(in) :: path !< NetCDF diagnostics output path.
    type(scaler_t), intent(in) :: scaler !< Fine-to-coarse grid scaler used for checkerboard diagnostics.
    integer(i4), intent(in) :: sub_map(:) !< Fine-river SCC label for each L0 node.
    logical, intent(in) :: leaving_cells(:) !< Fine-node mask of coarse-cell exits.
    logical, intent(in) :: stream_mask(:) !< Fine-node mask of retained coarse links.
    logical, intent(in) :: is_link_start(:) !< Fine-node mask of coarse-link starts.
    integer(i8), intent(in) :: lake_outlet_cells(:) !< Canonical L0 lake outlet nodes.
    integer(i8), intent(in) :: lake_ids(:) !< Stable IDs aligned with lake_outlet_cells.
    integer(i4), allocatable :: stream_sub(:)
    integer(i8), allocatable :: lake_outlet(:)
    logical, allocatable :: highlight(:)
    integer(i8) :: i

    call message("river_upscaler: write diagnostics to file: " // trim(path))
    allocate(stream_sub(this%fine_river%n_nodes))
    allocate(highlight(this%fine_river%n_nodes))
    allocate(lake_outlet(this%fine_river%n_nodes))
    !$omp parallel do default(shared)
    do i = 1_i8, this%fine_river%n_nodes
      stream_sub(i) = nodata_i4
      if (stream_mask(i)) stream_sub(i) = sub_map(i)
      highlight(i) = is_link_start(i) .or. this%fine_river%is_sink(i)
      lake_outlet(i) = 0_i8
    end do
    !$omp end parallel do
    if (size(lake_outlet_cells, kind=i8) /= size(lake_ids, kind=i8)) &
      call error_message("river_upscaler diagnostics: lake outlet node and ID counts differ")
    !$omp parallel do default(shared) schedule(static)
    do i = 1_i8, size(lake_ids, kind=i8)
      lake_outlet(lake_outlet_cells(i)) = lake_ids(i)
    end do
    !$omp end parallel do
    call this%fine_river%export( &
      path        = path, &
      sub_map     = sub_map, &
      leaving     = leaving_cells, &
      stream_mask = stream_mask, &
      stream_sub  = stream_sub, &
      lake_outlet = lake_outlet, &
      highlight   = highlight, &
      factor      = scaler%factor &
    )
  end subroutine river_upscaler_write_diagnostics

  !> \brief calculate the celerity c_i from slope s_i (i - cell index)
  subroutine river_upscaler_celerity(this, gamma, celerity, constant_celerity, slope)
    implicit none
    class(river_upscaler_t), target, intent(inout) :: this !< Initialized upscaler with retained fine stream mask.
    real(dp), intent(in) :: gamma !< model parameter: c_i = gamma * sqrt(s_i) or c = gamma
    real(dp), allocatable, intent(out) :: celerity(:) !< celerity of the link starting at each coarse node
    logical, optional, intent(in) :: constant_celerity !< whether celerity is assumed constant: c = gamma (default: .false.)
    real(dp), optional, intent(in) :: slope(:) !< [%] river slope on fine grid: size(fine\%ncells)
    integer(i8) :: i, cell
    real(dp) :: n
    real(dp), allocatable :: fine_celerity(:)

    if (optval(constant_celerity, .false.)) then
      call message("river_upscaler: constant celerity assumed, set c_i = gamma for all coarse nodes")
      ! constant celerity, no need to calculate from fine river
      call this%coarse_river%calc_celerity(gamma, celerity, constant_celerity)
      return
    end if

    if (.not.allocated(this%stream_mask)) then
      call error_message("river_upscaler%calc_celerity: variable celerity requires retain_stream_mask=.true. during init")
    end if

    ! first calculate celerity on fine river
    call this%fine_river%calc_celerity(gamma, fine_celerity, constant_celerity, slope, this%stream_mask)

    call message("river_upscaler: calculate celerity on coarse river from fine river")
    allocate(celerity(this%coarse_river%n_nodes))
    !$omp parallel do default(shared) private(i, cell, n)
    do i = 1_i8, this%coarse_river%n_nodes
      if (this%coarse_river%is_sink(i)) then
        celerity(i) = 1.0_dp
        cycle
      end if
      cell = this%link_start(i)
      celerity(i) = 0.0_dp
      n = 0.0_dp
      ! one pass algorithm for harmonic mean:
      ! 0. M  = 0               -> 0 as initial value for the [M]ean of inverses
      ! 1. M += (1/v - M) / n   -> updating the mean with the weighted deviation of new value (n as counter)
      ! 2. H  = 1 / M           -> [H]armonic mean is then the inverse of M
      do while (cell /= this%link_end(i))
        n = n + 1.0_dp
        celerity(i) = celerity(i) + ( 1.0_dp / fine_celerity(cell) - celerity(i) ) / n
        cell = this%fine_river%down(cell)
      end do
      ! finalize harmonic mean for celerity
      celerity(i) = 1.0_dp / celerity(i)
    end do
    !$omp end parallel do

  end subroutine river_upscaler_celerity

  subroutine river_upscaler_destroy(this)
    implicit none
    class(river_upscaler_t), target, intent(inout) :: this !< Upscaler whose temporary retained data is released.

    if (allocated(this%link_start)) deallocate(this%link_start)
    if (allocated(this%stream_mask)) deallocate(this%stream_mask)
    if (allocated(this%link_end)) deallocate(this%link_end)
    if (allocated(this%scc_coarse_gauges)) deallocate(this%scc_coarse_gauges)

  end subroutine river_upscaler_destroy

end module mo_river_upscaler
