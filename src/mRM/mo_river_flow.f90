!> \file    mo_river_flow.f90
!> \copydoc mo_river_flow

!> \brief   River discharge routing.
!> \details This module contains mRM-owned discharge state and a flow processor for
!! runoff preparation, lake exchange, and the shared river DAG traversal.
!> \version 0.1
!> \changelog
!! - 2026: discharge state and flow processor extracted from the river router.
!! - 2026: processor-owned scaling workspaces and cached lake topology.
!> \authors Sebastian Mueller
!> \date    October 2026
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_flow
  use mo_kind, only: i4, i8, dp
  use mo_river, only: river_t
  use mo_river_tools, only: river_muskingum_state_t, river_depth_to_flux
  use mo_grid, only: grid_t
  use mo_grid_scaler, only: scaler_t
  use mo_message, only: error_message
  use mo_utils, only: optval
  implicit none
  private
  public :: river_discharge_state_t, river_flow_t

  !> \class river_discharge_state_t
  !> \brief mRM-owned current and previous discharge fluxes [m3 s-1].
  type :: river_discharge_state_t
    real(dp), allocatable :: discharge(:) !< Current node inflow.
    real(dp), allocatable :: tributary(:) !< Current routed link outflow.
    real(dp), allocatable :: previous_discharge(:) !< Previous node inflow.
    real(dp), allocatable :: previous_tributary(:) !< Previous routed link outflow.
  contains
    procedure :: allocate => river_discharge_state_allocate
    procedure :: destroy => river_discharge_state_destroy
  end type river_discharge_state_t

  !> \class river_flow_t
  !> \brief Router-owned runoff and lake scratch used to route mRM-owned discharge state.
  type :: river_flow_t
    type(river_t), pointer :: river => null() !< Borrowed network; the mRM-owned target outlives this processor.
    integer(i4) :: n_lakes = 0_i4 !< Number of canonical lakes in exchange order.
    real(dp), allocatable :: scale_grid(:, :) !< Reusable unpacked runoff-volume workspace.
    logical :: lake_routing = .false. !< Whether canonical lake nodes are present.
    real(dp), allocatable :: acc_runoff(:) !< Accumulated L1 runoff [mm].
    real(dp), allocatable :: scaled_flux(:) !< L3 mapped runoff rate [m3 s-1].
    real(dp), allocatable :: scale_work(:) !< L1 volume work buffer.
    real(dp), allocatable :: runoff(:) !< Local node runoff [m3 s-1].
    integer(i4), allocatable :: lake_index(:) !< Canonical lake index for each node.
  contains
    procedure :: init => river_flow_init
    procedure :: destroy => river_flow_destroy
    procedure :: accumulate => river_flow_accumulate
    procedure :: prepare => river_flow_prepare
    procedure :: route_nodes => river_flow_route_nodes
    procedure :: advance => river_flow_advance
  end type river_flow_t
contains
  !> \brief Allocate and reset mRM-owned discharge histories.
  subroutine river_discharge_state_allocate(self, n_nodes)
    class(river_discharge_state_t), intent(inout), target :: self !< Discharge state to initialize.
    integer(i8), intent(in) :: n_nodes !< Number of routing nodes.
    call self%destroy()
    allocate(self%discharge(n_nodes), self%tributary(n_nodes), &
      self%previous_discharge(n_nodes), self%previous_tributary(n_nodes))
    self%discharge = 0.0_dp
    self%tributary = 0.0_dp
    self%previous_discharge = 0.0_dp
    self%previous_tributary = 0.0_dp
  end subroutine river_discharge_state_allocate

  !> \brief Release mRM-owned discharge histories.
  subroutine river_discharge_state_destroy(self)
    class(river_discharge_state_t), intent(inout), target :: self !< Discharge state to release.
    if (allocated(self%discharge)) deallocate(self%discharge)
    if (allocated(self%tributary)) deallocate(self%tributary)
    if (allocated(self%previous_discharge)) deallocate(self%previous_discharge)
    if (allocated(self%previous_tributary)) deallocate(self%previous_tributary)
  end subroutine river_discharge_state_destroy

  !> \brief Allocate runoff scratch and build the optional canonical lake map.
  subroutine river_flow_init(self, river, input_cells, runoff_cells, scc, needs_scale_work, lake_nodes, input_grid)
    class(river_flow_t), intent(inout), target :: self !< Flow processor to initialize.
    type(river_t), intent(in), target :: river !< Routing network.
    integer(i8), intent(in) :: input_cells !< L1 active cells.
    integer(i8), intent(in) :: runoff_cells !< L3 runoff support cells.
    logical, intent(in) :: scc !< Whether SCC topology is active.
    logical, intent(in) :: needs_scale_work !< Whether upscaling requires an L1 volume buffer.
    integer(i8), optional, intent(in) :: lake_nodes(:) !< Canonical lake nodes in point-set order.
    type(grid_t), optional, intent(in) :: input_grid !< L1 support for the unpacked scaling workspace.
    integer(i8) :: n

    call self%destroy()
    self%river => river
    allocate(self%acc_runoff(input_cells), self%runoff(river%n_nodes))
    if (scc) allocate(self%scaled_flux(runoff_cells))
    if (needs_scale_work) then
      allocate(self%scale_work(input_cells))
      if (present(input_grid)) allocate(self%scale_grid(input_grid%nx,input_grid%ny))
    end if
    self%acc_runoff = 0.0_dp
    self%runoff = 0.0_dp
    if (present(lake_nodes)) then
      if (size(lake_nodes, kind=i8) > 0_i8) then
        allocate(self%lake_index(river%n_nodes), source=0_i4)
        do n = 1_i8, size(lake_nodes, kind=i8)
          if (lake_nodes(n) < 1_i8 .or. lake_nodes(n) > river%n_nodes) &
            call error_message("river_flow: lake node is outside the river network")
          if (self%lake_index(lake_nodes(n)) /= 0_i4) call error_message("river_flow: duplicate canonical lake node")
          self%lake_index(lake_nodes(n)) = int(n, i4)
        end do
        self%n_lakes = int(size(lake_nodes), i4)
        self%lake_routing = .true.
      end if
    end if
  end subroutine river_flow_init

  !> \brief Release all router-owned flow scratch.
  subroutine river_flow_destroy(self)
    class(river_flow_t), intent(inout), target :: self !< Flow processor to release.
    if (allocated(self%acc_runoff)) deallocate(self%acc_runoff)
    if (allocated(self%scaled_flux)) deallocate(self%scaled_flux)
    if (allocated(self%scale_work)) deallocate(self%scale_work)
    if (allocated(self%runoff)) deallocate(self%runoff)
    if (allocated(self%lake_index)) deallocate(self%lake_index)
    if (allocated(self%scale_grid)) deallocate(self%scale_grid)
    self%n_lakes = 0_i4
    self%lake_routing = .false.
    nullify(self%river)
  end subroutine river_flow_destroy

  !> \brief Accumulate model-step runoff amounts until routing is due.
  subroutine river_flow_accumulate(self, input_runoff)
    class(river_flow_t), intent(inout), target :: self !< Flow processor with runoff accumulator.
    real(dp), intent(in) :: input_runoff(:) !< L1 runoff depth [mm].
    self%acc_runoff = self%acc_runoff + input_runoff
  end subroutine river_flow_accumulate

  !> \brief Map mean runoff to node fluxes and reset the completed-interval accumulator.
  subroutine river_flow_prepare(self, scaler, muskingum, input_step, node_runoff_cell)
    class(river_flow_t), intent(inout), target :: self !< Flow processor owning mapped runoff scratch.
    type(scaler_t), intent(inout) :: scaler !< Router-owned L1-to-L3 depth-to-flux map.
    type(river_muskingum_state_t), intent(in) :: muskingum !< Shared coefficients and completed routing cadence.
    integer(i4), intent(in) :: input_step !< Represented runoff support [h].
    integer(i8), intent(in), optional :: node_runoff_cell(:) !< Router-owned mapping from nodes to land runoff cells.
    integer(i8) :: n, c

    ! Average repeated runoff amounts across model updates before converting the represented input support to a rate.
    if (muskingum%accumulations > 1_i4) self%acc_runoff = self%acc_runoff / muskingum%accumulations
    ! distribute runoff in case of SCC
    if (self%river%scc) then
      call river_depth_to_flux(scaler, self%acc_runoff, self%scaled_flux, input_step, &
        self%scale_work, self%scale_grid)
      ! Keep SIMD but avoid the combined OpenMP 4.x form that NAG rejects here.
      !$omp parallel default(none) shared(self,node_runoff_cell) private(c, n)
      !$omp do simd schedule(static)
      do n = 1_i8, self%river%n_nodes
        if (present(node_runoff_cell)) then
          c = node_runoff_cell(n)
        else
          c = self%river%node_cell(n)
        end if
        if (c > 0_i8) then
          self%runoff(n) = self%scaled_flux(c) * self%river%area_fraction(n)
        else
          self%runoff(n) = 0.0_dp
        end if
      end do
      !$omp end do simd
      !$omp end parallel
    else
      call river_depth_to_flux(scaler, self%acc_runoff, self%runoff, input_step, self%scale_work, self%scale_grid)
    end if
    self%acc_runoff = 0.0_dp
  end subroutine river_flow_prepare

  !> \brief Route a complete level or consecutive serial levels in topological order.
  !> \details Parallel callers share a team and bounds; the loop barrier completes the range before returning.
  recursive subroutine river_flow_route_nodes(self, first, last, state, muskingum, lake_outflow, lake_inflow, parallel)
    class(river_flow_t), intent(in), target :: self !< Prepared local runoff and lake map.
    integer(i8), intent(in) :: first !< First topological level to route.
    integer(i8), intent(in) :: last !< Last topological level to route.
    type(river_discharge_state_t), intent(inout), target :: state !< Discharge fluxes to update.
    type(river_muskingum_state_t), intent(in) :: muskingum !< Read-only routing coefficients and cadence.
    real(dp), optional, intent(in) :: lake_outflow(:) !< Prescribed lake release [m3 s-1].
    real(dp), optional, intent(inout) :: lake_inflow(:) !< Collected lake inflow [m3 s-1].
    logical, optional, intent(in) :: parallel !< Collective static worksharing; default is a serial topological traversal.
    integer(i8) :: j, first_node, last_node

    first_node = self%river%order%level_start(first)
    last_node = self%river%order%level_end(last)

    if (optval(parallel,.false.)) then
      ! Complete discharge throughout this range before dependent routing processors start.
      !$omp do schedule(static)
      do j = first_node, last_node
        call river_flow_process_node(j)
      end do
      !$omp end do
    else
      do j = first_node, last_node
        call river_flow_process_node(j)
      end do
    end if
  contains
    !> \brief Apply the node balance using the enclosing range's shared state and coefficient arrays.
    subroutine river_flow_process_node(position)
      integer(i8), intent(in) :: position !< Position in the node order.
      integer(i8) :: n
      integer(i8), pointer, contiguous :: edges(:)
      integer(i4) :: lake, m
      real(dp) :: inflow, prev_inflow, prev_routed, nu1, nu2

      n = self%river%order%id(position)
      inflow = self%runoff(n)
      call self%river%src_view(n, edges)
      !$omp simd reduction(+:inflow)
      do m = 1_i4, size(edges)
        inflow = inflow + state%tributary(edges(m))
      end do
      lake = 0_i4
      if (allocated(self%lake_index)) lake = self%lake_index(n)
      if (lake > 0_i4) then
        if (.not.present(lake_outflow) .or. .not.present(lake_inflow)) &
          call error_message("river_flow: lake routing requires inflow and outflow arrays")
        lake_inflow(lake) = lake_inflow(lake) + inflow
        inflow = lake_outflow(lake)
      end if
      state%discharge(n) = inflow
      if (self%river%is_sink(n)) then
        state%tributary(n) = 0.0_dp
        return
      end if
      if (allocated(self%lake_index)) then
        if (lake == 0_i4 .and. self%lake_index(self%river%down(n)) > 0_i4) then
          state%tributary(n) = inflow
          return
        end if
      end if
      nu1 = muskingum%nu1(n)
      nu2 = muskingum%nu2(n)
      prev_routed = state%previous_tributary(n)
      prev_inflow = state%previous_discharge(n)
      state%tributary(n) = prev_routed + nu1 * (prev_inflow - prev_routed) + nu2 * (inflow - prev_inflow)
    end subroutine river_flow_process_node
  end subroutine river_flow_route_nodes

  !> \brief Swap discharge history buffers without copying them.
  subroutine river_flow_advance(self, state)
    class(river_flow_t), intent(in), target :: self !< Flow processor.
    type(river_discharge_state_t), intent(inout), target :: state !< Discharge fluxes to advance.
    real(dp), allocatable :: tmp(:)
    call move_alloc(state%previous_discharge, tmp)
    call move_alloc(state%discharge, state%previous_discharge)
    call move_alloc(tmp, state%discharge)
    call move_alloc(state%previous_tributary, tmp)
    call move_alloc(state%tributary, state%previous_tributary)
    call move_alloc(tmp, state%tributary)
  end subroutine river_flow_advance
end module mo_river_flow
