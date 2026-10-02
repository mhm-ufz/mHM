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
  use mo_constants, only: nodata_dp
  use mo_netcdf, only: NcDataset, NcVariable, NcDimension
  use mo_river, only: river_t
  use mo_river_tools, only: river_muskingum_state_t, river_depth_to_flux
  use mo_grid, only: grid_t
  use mo_grid_scaler, only: scaler_t
  use mo_message, only: error_message, message
  use mo_utils, only: optval
  implicit none
  private
  public :: river_flow_state_t, river_flow_t

  !> \class river_flow_state_t
  !> \brief mRM-owned flow histories and stable published means [m3 s-1].
  type :: river_flow_state_t
    real(dp), allocatable :: previous_input(:) !< Previous node inflow.
    real(dp), allocatable :: previous_routed(:) !< Previous routed link outflow.
    real(dp), allocatable :: current_input(:) !< Current node inflow.
    real(dp), allocatable :: current_routed(:) !< Current routed link outflow.
    real(dp), allocatable :: discharge(:) !< Completed mean node discharge.
    real(dp), allocatable :: lake_inflow(:) !< Completed mean inflow in canonical lake order.
  contains
    procedure :: allocate => river_flow_state_allocate
    procedure :: reset => river_flow_state_reset
    procedure :: destroy => river_flow_state_destroy
  end type river_flow_state_t

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
    real(dp), allocatable :: sum_discharge(:) !< Sum of node discharge over numerical substeps [m3 s-1].
    real(dp), allocatable :: sum_lake_inflow(:) !< Sum of canonical lake inflow over numerical substeps [m3 s-1].
    integer(i4), allocatable :: lake_index(:) !< Canonical lake index for each node.
  contains
    procedure :: init => river_flow_init
    procedure :: to_restart_dataset => river_flow_to_restart_dataset
    procedure :: destroy => river_flow_destroy
    procedure :: accumulate => river_flow_accumulate
    procedure :: prepare => river_flow_prepare
    procedure :: route_nodes => river_flow_route_nodes
    procedure :: advance => river_flow_advance
    procedure :: complete => river_flow_complete
  end type river_flow_t
contains
  !> \brief Allocate mRM-owned histories and published flow arrays before exchange binding.
  subroutine river_flow_state_allocate(self, n_nodes, n_lakes)
    class(river_flow_state_t), intent(inout), target :: self !< Flow state to initialize.
    integer(i8), intent(in) :: n_nodes !< Number of routing nodes.
    integer(i4), optional, intent(in) :: n_lakes !< Number of canonical lakes; zero or absent disables lake publication.
    call self%destroy()
    allocate(self%previous_input(n_nodes), self%previous_routed(n_nodes), &
      self%current_input(n_nodes), self%current_routed(n_nodes), self%discharge(n_nodes))
    if (present(n_lakes)) then
      if (n_lakes > 0_i4) allocate(self%lake_inflow(n_lakes))
    end if
    call self%reset()
  end subroutine river_flow_state_allocate

  !> \brief Reset histories and completed means without replacing exchange targets.
  subroutine river_flow_state_reset(self)
    class(river_flow_state_t), intent(inout), target :: self !< Flow state to reset in place.
    if (allocated(self%current_input)) self%current_input = 0.0_dp
    if (allocated(self%current_routed)) self%current_routed = 0.0_dp
    if (allocated(self%previous_input)) self%previous_input = 0.0_dp
    if (allocated(self%previous_routed)) self%previous_routed = 0.0_dp
    if (allocated(self%discharge)) self%discharge = 0.0_dp
    if (allocated(self%lake_inflow)) self%lake_inflow = 0.0_dp
  end subroutine river_flow_state_reset

  !> \brief Release mRM-owned histories and published flow arrays after exchange detachment.
  subroutine river_flow_state_destroy(self)
    class(river_flow_state_t), intent(inout), target :: self !< Flow state to release.
    if (allocated(self%current_input)) deallocate(self%current_input)
    if (allocated(self%current_routed)) deallocate(self%current_routed)
    if (allocated(self%previous_input)) deallocate(self%previous_input)
    if (allocated(self%previous_routed)) deallocate(self%previous_routed)
    if (allocated(self%discharge)) deallocate(self%discharge)
    if (allocated(self%lake_inflow)) deallocate(self%lake_inflow)
  end subroutine river_flow_state_destroy

  !> \brief Initialize flow histories and scratch, optionally restoring previous fluxes from restart.
  subroutine river_flow_init(self, state, river, input_cells, runoff_cells, scc, needs_scale_work, lake_nodes, input_grid, nc, read_fluxes)
    class(river_flow_t), intent(inout), target :: self !< Flow processor to initialize.
    type(river_flow_state_t), intent(inout), target :: state !< mRM-owned flow state.
    type(river_t), intent(in), target :: river !< Routing network.
    integer(i8), intent(in) :: input_cells !< L1 active cells.
    integer(i8), intent(in) :: runoff_cells !< L3 runoff support cells.
    logical, intent(in) :: scc !< Whether SCC topology is active.
    logical, intent(in) :: needs_scale_work !< Whether upscaling requires an L1 volume buffer.
    integer(i8), optional, intent(in) :: lake_nodes(:) !< Canonical lake nodes in point-set order.
    type(grid_t), optional, intent(in) :: input_grid !< L1 support for the unpacked scaling workspace.
    type(NcDataset), intent(in), optional :: nc !< Open restart dataset replacing fresh flow histories.
    logical, intent(in), optional :: read_fluxes !< Restore previous fluxes; false resets them (default true).
    type(NcVariable) :: var
    integer(i8) :: n

    call self%destroy()
    self%river => river
    if (.not.allocated(state%current_input)) call error_message("river_flow: mRM-owned flow state is not allocated")
    if (size(state%current_input, kind=i8) /= river%n_nodes) &
      call error_message("river_flow: flow state does not match the river network")
    ! Current buffers are overwritten on the next substep; histories start cold unless explicitly restored.
    state%current_input = 0.0_dp
    state%current_routed = 0.0_dp
    state%previous_input = 0.0_dp
    state%previous_routed = 0.0_dp
    if (present(nc)) then
      if (optval(read_fluxes, .true.)) then
        ! Preserve compatibility with older restarts that omitted either history.
        if (nc%hasVariable("previous_discharge")) then
          var = nc%getVariable("previous_discharge")
          call var%readInto(state%previous_input)
        end if
        if (nc%hasVariable("previous_tributary")) then
          var = nc%getVariable("previous_tributary")
          call var%readInto(state%previous_routed)
        end if
      end if
    end if
    allocate(self%acc_runoff(input_cells), self%runoff(river%n_nodes), self%sum_discharge(river%n_nodes))
    if (scc) allocate(self%scaled_flux(runoff_cells))
    if (needs_scale_work) then
      allocate(self%scale_work(input_cells))
      if (present(input_grid)) allocate(self%scale_grid(input_grid%nx,input_grid%ny))
    end if
    self%acc_runoff = 0.0_dp
    self%runoff = 0.0_dp
    self%sum_discharge = 0.0_dp
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
        allocate(self%sum_lake_inflow(self%n_lakes), source=0.0_dp)
      end if
    end if
  end subroutine river_flow_init

  !> \brief Persist flow histories using the existing restart variable names and metadata.
  subroutine river_flow_to_restart_dataset(self, state, nc, node_dim)
    class(river_flow_t), intent(in), target :: self !< Flow processor.
    type(river_flow_state_t), intent(in), target :: state !< mRM-owned flow histories.
    type(NcDataset), intent(inout) :: nc !< Open restart dataset.
    type(NcDimension), intent(in) :: node_dim !< River node dimension.
    type(NcVariable) :: nc_var

    ! Write only the previous arrays: the current arrays are scratch buffers that
    ! are overwritten before they are read by the next routing calculation.
    if ( allocated(state%previous_input) ) then
      call message("writing previous_discharge to restart_file")
      nc_var = nc%setVariable("previous_discharge", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "previous_discharge")
      call nc_var%setAttribute("units", "m3 s-1")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(state%previous_input)
    end if

    ! write previous_tributary
    if ( allocated(state%previous_routed) ) then
      call message("writing previous_tributary to restart_file")
      nc_var = nc%setVariable("previous_tributary", "f64", [node_dim])
      call nc_var%setAttribute("long_name", "previous_tributary")
      call nc_var%setAttribute("units", "m3 s-1")
      call nc_var%setFillValue(nodata_dp)
      call nc_var%setAttribute("missing_value", nodata_dp)
      call nc_var%setData(state%previous_routed)
    end if

  end subroutine river_flow_to_restart_dataset

  !> \brief Release all router-owned flow scratch.
  subroutine river_flow_destroy(self)
    class(river_flow_t), intent(inout), target :: self !< Flow processor to release.
    if (allocated(self%acc_runoff)) deallocate(self%acc_runoff)
    if (allocated(self%scaled_flux)) deallocate(self%scaled_flux)
    if (allocated(self%scale_work)) deallocate(self%scale_work)
    if (allocated(self%runoff)) deallocate(self%runoff)
    if (allocated(self%lake_index)) deallocate(self%lake_index)
    if (allocated(self%scale_grid)) deallocate(self%scale_grid)
    if (allocated(self%sum_discharge)) deallocate(self%sum_discharge)
    if (allocated(self%sum_lake_inflow)) deallocate(self%sum_lake_inflow)
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
    self%sum_discharge = 0.0_dp
    if (allocated(self%sum_lake_inflow)) self%sum_lake_inflow = 0.0_dp
  end subroutine river_flow_prepare

  !> \brief Route a complete level or consecutive serial levels in topological order.
  !> \details Parallel callers share a team and bounds; the loop barrier completes the range before returning.
  recursive subroutine river_flow_route_nodes(self, first, last, state, muskingum, lake_outflow, parallel)
    class(river_flow_t), intent(inout), target :: self !< Prepared local runoff and lake map.
    integer(i8), intent(in) :: first !< First topological level to route.
    integer(i8), intent(in) :: last !< Last topological level to route.
    type(river_flow_state_t), intent(inout), target :: state !< Discharge fluxes to update.
    type(river_muskingum_state_t), intent(in) :: muskingum !< Read-only routing coefficients and cadence.
    real(dp), optional, intent(in) :: lake_outflow(:) !< Prescribed lake release [m3 s-1].
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
        inflow = inflow + state%current_routed(edges(m))
      end do
      lake = 0_i4
      if (allocated(self%lake_index)) lake = self%lake_index(n)
      if (lake > 0_i4) then
        if (.not.present(lake_outflow)) call error_message("river_flow: lake routing requires prescribed outflow")
        self%sum_lake_inflow(lake) = self%sum_lake_inflow(lake) + inflow
        inflow = lake_outflow(lake)
      end if
      state%current_input(n) = inflow
      self%sum_discharge(n) = self%sum_discharge(n) + inflow
      if (self%river%is_sink(n)) then
        state%current_routed(n) = 0.0_dp
        return
      end if
      if (allocated(self%lake_index)) then
        if (lake == 0_i4 .and. self%lake_index(self%river%down(n)) > 0_i4) then
          state%current_routed(n) = inflow
          return
        end if
      end if
      nu1 = muskingum%nu1(n)
      nu2 = muskingum%nu2(n)
      prev_routed = state%previous_routed(n)
      prev_inflow = state%previous_input(n)
      state%current_routed(n) = prev_routed + nu1 * (prev_inflow - prev_routed) + nu2 * (inflow - prev_inflow)
    end subroutine river_flow_process_node
  end subroutine river_flow_route_nodes

  !> \brief Swap discharge history buffers without copying them.
  subroutine river_flow_advance(self, state)
    class(river_flow_t), intent(in), target :: self !< Flow processor.
    type(river_flow_state_t), intent(inout), target :: state !< Discharge fluxes to advance.
    real(dp), allocatable :: tmp(:)
    call move_alloc(state%previous_input, tmp)
    call move_alloc(state%current_input, state%previous_input)
    call move_alloc(tmp, state%current_input)
    call move_alloc(state%previous_routed, tmp)
    call move_alloc(state%current_routed, state%previous_routed)
    call move_alloc(tmp, state%current_routed)
  end subroutine river_flow_advance

  !> \brief Publish completed discharge and canonical lake-inflow means without moving their targets.
  subroutine river_flow_complete(self, state, muskingum)
    class(river_flow_t), intent(in), target :: self !< Flow processor holding completed substep sums.
    type(river_flow_state_t), intent(inout), target :: state !< Flow output state.
    type(river_muskingum_state_t), intent(in) :: muskingum !< Number of substeps in the completed interval.
    state%discharge = self%sum_discharge / muskingum%iterations
    if (self%lake_routing) state%lake_inflow = self%sum_lake_inflow / muskingum%iterations
  end subroutine river_flow_complete
end module mo_river_flow
