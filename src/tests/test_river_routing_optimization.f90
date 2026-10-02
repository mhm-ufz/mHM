!> \file test_river_routing_optimization.f90
!> \brief Check forcing caches, remapping, dry/rewetting routing, and OpenMP level eligibility.
program test_river_routing_optimization
  use mo_domain, only: domain_t
  use mo_kind, only: i2, i4, i8, dp
  use mo_constants, only: T0_dp, nodata_dp
  use mo_grid, only: grid_t, bottom_up, top_down
  use mo_grid_scaler, only: scaler_t, up_a_mean, up_sum, down_nearest, down_scaling
  use mo_river, only: river_t, sink, d8_E
  use mo_river_router, only: river_router_t
  use mo_river_flow, only: river_flow_state_t, river_flow_t
  use mo_river_tools, only: river_muskingum_state_t
  use mo_river_temperature, only: river_thermal_t, river_thermal_state_t, water_heat_capacity_volume
  !$ use omp_lib, only: omp_set_num_threads, omp_set_dynamic, omp_get_max_threads
  implicit none
  integer(i4) :: mode
  !$ call omp_set_dynamic(.false.)
  !$ call omp_set_num_threads(2)
  do mode = 1_i4, 3_i4
    call check_mapping(mode)
  end do
  call check_levels(.false.)
  call check_levels(.true.)
  call check_level_worksharing()
  call check_accumulation_cadence()
  call check_flow_completion()
  call check_shared_node_map()
  call check_sealed_fraction()
contains
  !> \brief Compare optimized forcing/flux mapping against the existing scaler, including orientation and masks.
  subroutine check_mapping(mode)
    integer(i4), intent(in) :: mode !< Downscaling, equal grid, or upscaling fixture.
    type(grid_t), target :: fine, coarse, equal
    type(grid_t), pointer :: input_grid, river_grid
    type(river_t), target :: river
    type(river_t), pointer :: river_ptr
    type(river_router_t) :: router
    type(river_flow_state_t) :: flow
    type(river_thermal_state_t) :: state
    type(scaler_t) :: reference
    integer(i2), allocatable :: fdir(:)
    real(dp), allocatable :: input(:), air(:), shortwave(:), longwave(:), depth(:), expected(:)
    real(dp) :: parameters(5), expected_ta, expected_r0
    integer(i8) :: c, n, downstream, ci, co
    logical :: mask(4,2), completed

    mask = .true.
    mask(2,1) = .false.
    call fine%init(4_i4,2_i4,cellsize=2.0_dp,mask=mask,y_direction=bottom_up)
    call equal%init(4_i4,2_i4,cellsize=2.0_dp,mask=mask(:,2:1:-1),y_direction=top_down)
    call coarse%init(2_i4,1_i4,cellsize=4.0_dp,y_direction=top_down)
    ! Nonuniform physical areas exercise normalized area weights, independently of grid cellsize.
    do c = 1_i8, fine%ncells
      fine%cell_area(c) = 2.0_dp + real(c,dp)
    end do
    select case(mode)
    case(1)
      input_grid => coarse
      river_grid => fine
    case(2)
      input_grid => fine
      river_grid => equal
    case(3)
      input_grid => fine
      river_grid => coarse
    end select
    allocate(fdir(river_grid%ncells),source=d8_E)
    do c = 1_i8, river_grid%ncells
      if (river_grid%cell_ij(c,1) == river_grid%nx) fdir(c) = sink
    end do
    call river%from_fdir(fdir,river_grid)
    river_ptr => river
    call flow%allocate(river%n_nodes)
    call state%allocate(river%n_nodes)
    parameters = [10.0_dp,0.08_dp,1.26_dp,0.98_dp,20.0_dp]
    state%active = .true.
    state%parameters = parameters
    call router%init(river_ptr,flow,state,spread(1.0_dp,1,int(river%n_nodes)),input_grid, &
      input_step=1_i4,max_route_step=3600.0_dp,omp_level_thresh=0_i8)
    allocate(input(input_grid%ncells),air(input_grid%ncells),shortwave(input_grid%ncells), &
      longwave(input_grid%ncells),depth(input_grid%ncells),expected(river_grid%ncells))
    do c = 1_i8, input_grid%ncells
      air(c) = 5.0_dp + real(c,dp)
      shortwave(c) = 30.0_dp*real(c,dp)
      longwave(c) = 250.0_dp + real(c,dp)
    end do
    call reference%init(input_grid,river_grid,upscaling_operator=up_a_mean,downscaling_operator=down_nearest)
    call router%thermal%accumulate(air,air,shortwave,longwave)
    call router%thermal%prepare(router%scaler,router%muskingum,router%input_step,router%node_runoff_cell)
    call reference%execute(air,expected)
    call assert_true(all(abs(expected-router%thermal%air)<1.e-12_dp),'mapped air temperature')
    call reference%execute(shortwave,expected)
    call assert_true(all(abs(expected-router%thermal%ssrd)<1.e-12_dp),'mapped shortwave')
    call reference%execute(longwave,expected)
    call assert_true(all(abs(expected-router%thermal%strd)<1.e-12_dp),'mapped longwave')
    do n = 1_i8, river%n_nodes
      downstream = n
      if (.not.river%is_sink(n)) downstream = river%down(n)
      ci = n
      co = downstream
      if (allocated(river%node_cell)) then
        ci = river%node_cell(n)
        co = river%node_cell(downstream)
      end if
      expected_ta = 0.5_dp*(router%thermal%air(ci)+router%thermal%air(co))+T0_dp
      expected_r0 = (1.0_dp-parameters(2))*0.5_dp*(router%thermal%ssrd(ci)+router%thermal%ssrd(co)) + &
        0.5_dp*(router%thermal%strd(ci)+router%thermal%strd(co))
      call assert_true(abs(expected_ta-router%thermal%forcing%ta(n))<1.e-12_dp,'endpoint temperature cache')
      call assert_true(abs(expected_r0-router%thermal%forcing%r0(n))<1.e-12_dp,'endpoint radiation cache')
    end do
    ! A cache must refresh when meteorology changes in the following routing interval.
    air = air+3.0_dp
    call router%thermal%accumulate(air,air,shortwave,longwave)
    call router%thermal%prepare(router%scaler,router%muskingum,router%input_step,router%node_runoff_cell)
    call reference%execute(air,expected)
    call assert_true(all(abs(expected-router%thermal%air)<1.e-12_dp),'forcing refresh')
    input = 0.0_dp
    depth = 0.0_dp
    call router%update(input,flow,state,completed,temperature_depth=depth, &
      air_temperature=air,shortwave=shortwave,longwave=longwave)
    call assert_true(completed,'completed dry interval')
    call assert_true(all(state%temperature==nodata_dp) .and. all(state%energy_flux==0.0_dp),'dry publication')
    input = 0.1_dp
    depth = input*280.0_dp
    call router%update(input,flow,state,completed,temperature_depth=depth, &
      air_temperature=air,shortwave=shortwave,longwave=longwave)
    call assert_true(all(abs(state%temperature-(280.0_dp-T0_dp))<1.e-6_dp),'rewetting preserves energy temperature')
    call reference%init(input_grid,river_grid,upscaling_operator=up_sum,downscaling_operator=down_nearest)
    if (reference%scaling_mode == down_scaling) then
      call reference%downscale_nearest(input,expected)
      expected = expected*river_grid%cell_area
    else
      call reference%execute(input*input_grid%cell_area,expected)
    end if
    expected = expected*(0.001_dp/3600.0_dp)
    call assert_true(all(expected==router%flow%runoff),'depth-to-flux mapping is unchanged')
    call assert_true(all(abs(router%thermal%local_energy-water_heat_capacity_volume*expected*280.0_dp) < &
      1.e-10_dp*max(1.0_dp,maxval(abs(router%thermal%local_energy)))),'mapped lateral energy')
    ! Average repeated amounts before applying the shared depth-to-flux map.
    input = 0.1_dp
    depth = input*278.0_dp
    call router%flow%accumulate(input)
    call router%thermal%accumulate(depth,air,shortwave,longwave)
    input = 0.3_dp
    depth = input*282.0_dp
    call router%flow%accumulate(input)
    call router%thermal%accumulate(depth,air,shortwave,longwave)
    call assert_true(all(abs(router%flow%acc_runoff-0.4_dp)<1.e-12_dp),'both model updates accumulated')
    router%muskingum%accumulations = 2_i4
    call router%flow%prepare(router%scaler,router%muskingum,router%input_step,router%node_runoff_cell)
    call router%thermal%prepare(router%scaler,router%muskingum,router%input_step,router%node_runoff_cell)
    input = 0.2_dp
    if (reference%scaling_mode == down_scaling) then
      call reference%downscale_nearest(input,expected)
      expected = expected*river_grid%cell_area
    else
      call reference%execute(input*input_grid%cell_area,expected)
    end if
    expected = expected*(0.001_dp/3600.0_dp)
    call assert_true(all(abs(expected-router%flow%runoff)<1.e-12_dp),'averaged flow preparation')
    call assert_true(all(abs(router%thermal%local_energy-water_heat_capacity_volume*expected*281.0_dp) < &
      1.e-10_dp*max(1.0_dp,maxval(abs(router%thermal%local_energy)))),'averaged thermal preparation')
    call assert_true(all(router%flow%acc_runoff==0.0_dp) .and. all(router%thermal%acc_depth==0.0_dp),'preparation resets accumulation')
    input = 0.4_dp
    call router%flow%accumulate(input)
    call assert_true(all(router%flow%acc_runoff==input),'new interval starts from the reset runoff accumulator')

    ! A supplied inactive state must never activate or retain thermal scratch.
    state%active = .false.
    call router%init(river_ptr,flow,state,spread(1.0_dp,1,int(river%n_nodes)),input_grid, &
      input_step=1_i4,max_route_step=3600.0_dp,omp_level_thresh=0_i8)
    call assert_true(.not.state%active .and. .not.associated(router%thermal%river),'inactive thermal state stays inactive')
    call assert_true(.not.allocated(router%thermal%acc_depth),'inactive initialization releases thermal scratch')
    call router%destroy()
    call assert_true(allocated(flow%discharge) .and. allocated(state%temperature),'processor destruction preserves state')
    call assert_true(.not.allocated(router%thermal%forcing%ta),'forcing cache destruction')
    call flow%destroy()
    call state%destroy()
    call reference%destroy()
    call river%destroy()
  end subroutine check_mapping

  !> \brief Check exact threshold boundaries, root ordering, single levels, and one-thread routing.
  subroutine check_levels(all_sinks)
    logical, intent(in) :: all_sinks !< Use a single-level sink network instead of horizontal chains.
    type(grid_t), target :: grid
    type(grid_t), pointer :: grid_ptr
    type(river_t), target :: river
    type(river_t), pointer :: river_ptr
    type(river_router_t) :: router
    type(river_flow_state_t) :: flow
    type(river_thermal_state_t) :: thermal
    integer(i2) :: fdir(32)
    integer(i8) :: threshold, c
    integer(i4) :: threads
    call grid%init(8_i4,4_i4,cellsize=2000.0_dp)
    grid_ptr => grid
    fdir = d8_E
    do c=1_i8,grid%ncells
      if (grid%cell_ij(c,1)==grid%nx .or. all_sinks) fdir(c)=sink
    end do
    call river%from_fdir(fdir,grid_ptr)
    river_ptr => river
    call flow%allocate(river%n_nodes)
    threads = 1_i4
    !$ threads = omp_get_max_threads()
    do threshold = 0_i8, maxval([4_i8,river%n_nodes+1_i8])
      call router%init(river_ptr,flow,thermal,spread(1.0_dp,1,32),grid_ptr, &
        omp_level_thresh=threshold,root_levels=.true.)
      if (threshold == 0_i8 .or. threads == 1_i4) then
        call assert_true(.not.any(router%parallel_level),'serial level selection')
      else
        call assert_true(all(router%parallel_level .eqv. (river%order%level_size>=threshold)),'threshold level eligibility')
      end if
    end do
    !$ call omp_set_num_threads(1)
    call router%init(river_ptr,flow,thermal,spread(1.0_dp,1,32),grid_ptr,omp_level_thresh=1_i8)
    call assert_true(.not.any(router%parallel_level),'one thread disables node-pass parallelism')
    !$ call omp_set_num_threads(2)
    call router%destroy()
    call flow%destroy()
    call river%destroy()
  end subroutine check_levels

  !> \brief Compare workshared routing through confluences, serial gaps, and tails with a serial reference.
  subroutine check_level_worksharing()
    type(grid_t), target :: grid
    type(grid_t), pointer :: grid_ptr
    type(river_t), target :: river
    type(river_t), pointer :: river_ptr
    type(river_router_t) :: reference, workshared
    type(river_flow_state_t) :: reference_flow, workshared_flow
    type(river_thermal_state_t) :: reference_thermal, workshared_thermal
    integer(i8) :: down(16), threshold, n
    integer(i4) :: thermal_mode, ordering, step, threads
    real(dp) :: input(16), depth(16), air(16), sw(16), lw(16)
    real(dp) :: parameters(5)
    logical :: reference_completed, workshared_completed

    ! A ten-node main stem has tributaries at two depths, producing separated wide root-order levels.
    down = [2_i8,3_i8,4_i8,5_i8,6_i8,7_i8,8_i8,9_i8,10_i8,0_i8,3_i8,3_i8,4_i8,7_i8,7_i8,8_i8]
    call grid%init(16_i4,1_i4,cellsize=2000.0_dp)
    grid_ptr => grid
    call river%init(down)
    river%grid => grid
    allocate(river%is_sink(16), source=down==0_i8)
    allocate(river%link_length(16), source=1000.0_dp)
    where (river%is_sink) river%link_length=0.0_dp
    river_ptr => river
    parameters = [10.0_dp,0.08_dp,1.26_dp,0.98_dp,20.0_dp]
    threads = 1_i4
    !$ threads = omp_get_max_threads()
    do ordering = 0_i4, 1_i4
      call river%calc_order(root=ordering==1_i4,sort=.true.)
      do thermal_mode = 0_i4, 1_i4
        do threshold = 1_i8, 3_i8, 2_i8
          call reference_flow%allocate(river%n_nodes)
          call workshared_flow%allocate(river%n_nodes)
          if (thermal_mode==1_i4) then
            call reference_thermal%allocate(river%n_nodes)
            call workshared_thermal%allocate(river%n_nodes)
            reference_thermal%active = .true.
            reference_thermal%parameters = parameters
            reference_thermal%area = 25.0_dp
            workshared_thermal%active = .true.
            workshared_thermal%parameters = parameters
            workshared_thermal%area = 25.0_dp
          end if
          call reference%init(river_ptr,reference_flow,reference_thermal,spread(1.0_dp,1,16),grid_ptr,omp_level_thresh=0_i8)
          call workshared%init(river_ptr,workshared_flow,workshared_thermal,spread(1.0_dp,1,16),grid_ptr,omp_level_thresh=threshold)
          if (ordering==1_i4 .and. threshold==3_i8 .and. threads>1_i4) then
            call assert_true(.not.workshared%parallel_level(1) .and. workshared%parallel_level(6) .and. .not.any(workshared%parallel_level(7:)),'serial prefix, last eligible level and tail')
            call assert_true(all(workshared%serial_level_end(3:5)==5_i8),'consecutive intermediate serial levels share one section')
            call assert_true(all(workshared%serial_level_end(7:10)==10_i8),'serial tail spans all remaining levels')
          end if
          do step = 1_i4, 8_i4
            do n = 1_i8, 16_i8
              input(n) = 0.1_dp+0.005_dp*real(n,dp)+0.001_dp*real(step,dp)
              depth(n) = input(n)*(277.0_dp+0.1_dp*real(n,dp))
              air(n) = 5.0_dp+0.2_dp*real(n,dp)
              sw(n) = 150.0_dp*real(mod(step,2_i4),dp)
              lw(n) = 300.0_dp
            end do
            call reference%update(input,reference_flow,reference_thermal,reference_completed, &
              temperature_depth=depth,air_temperature=air,shortwave=sw,longwave=lw)
            call workshared%update(input,workshared_flow,workshared_thermal,workshared_completed, &
              temperature_depth=depth,air_temperature=air,shortwave=sw,longwave=lw)
            call assert_true(reference_completed .eqv. workshared_completed,'workshared completion cadence')
            call assert_true(all(abs(reference_flow%discharge-workshared_flow%discharge)<1.e-12_dp),'workshared confluence discharge')
            call assert_true(all(abs(reference_flow%previous_routed-workshared_flow%previous_routed)<1.e-12_dp),'workshared Muskingum histories')
            if (thermal_mode==1_i4) then
              call assert_true(all(abs(reference_thermal%temperature-workshared_thermal%temperature)<1.e-6_dp),'flow before temperature within every level')
              call assert_true(all(abs(reference_thermal%energy_flux-workshared_thermal%energy_flux)<1.e-3_dp+1.e-8_dp*abs(reference_thermal%energy_flux)),'workshared thermal energy conservation')
            end if
          end do
          call reference%destroy()
          call workshared%destroy()
          call reference_flow%destroy()
          call workshared_flow%destroy()
          call reference_thermal%destroy()
          call workshared_thermal%destroy()
        end do
      end do
    end do
    call river%destroy()
  end subroutine check_level_worksharing

  !> \brief Check shared cadence, publication, and energy conservation across two model updates.
  subroutine check_accumulation_cadence()
    type(grid_t), target :: grid
    type(grid_t), pointer :: grid_ptr
    type(river_t), target :: river
    type(river_t), pointer :: river_ptr
    type(river_router_t) :: router
    type(river_flow_state_t) :: flow
    type(river_thermal_state_t) :: thermal
    real(dp) :: input(2), depth(2), air(2), sw(2), lw(2), previous_discharge(2)
    integer(i4) :: step
    logical :: completed

    call grid%init(2_i4,1_i4,cellsize=16000.0_dp)
    grid_ptr => grid
    call river%from_fdir([d8_E,sink],grid_ptr)
    river_ptr => river
    call flow%allocate(river%n_nodes)
    call thermal%allocate(river%n_nodes)
    thermal%active = .true.
    thermal%parameters = [10.0_dp,0.08_dp,1.26_dp,0.98_dp,20.0_dp]
    call router%init(river_ptr,flow,thermal,[1.0_dp,1.0_dp],grid_ptr, &
      input_step=1_i4,model_step=1_i4,max_route_step=7200.0_dp,omp_level_thresh=0_i8)
    call assert_true(router%muskingum%accumulations==2_i4,'two model updates per routing interval')
    air = 5.0_dp
    sw = 0.0_dp
    lw = 250.0_dp
    do step = 1_i4, 4_i4
      previous_discharge = flow%discharge
      input = 0.1_dp*step
      depth = input*280.0_dp
      call router%update(input,flow,thermal,completed, &
        temperature_depth=depth,air_temperature=air,shortwave=sw,longwave=lw)
      if (mod(step,2_i4)==1_i4) then
        call assert_true(.not.completed .and. router%input_count==1_i4,'router holds a partial interval')
        call assert_true(all(flow%discharge==previous_discharge),'partial interval preserves published discharge')
      else
        call assert_true(completed .and. router%input_count==0_i4,'router completes and resets its interval')
        call assert_true(all(router%flow%acc_runoff==0.0_dp) .and. all(router%thermal%acc_depth==0.0_dp),'both processors reset accumulators')
        call assert_true(all(abs(thermal%temperature-(280.0_dp-T0_dp))<1.e-6_dp),'two-update thermal mixing')
        call assert_true(all(abs(thermal%energy_flux-water_heat_capacity_volume*flow%discharge*280.0_dp) < &
          1.e-3_dp+1.e-8_dp*abs(thermal%energy_flux)),'two-update thermal energy conservation')
      end if
    end do
    call router%destroy()
    call assert_true(router%input_count==0_i4 .and. router%input_step==0_i4,'destroy resets shared cadence')
    call flow%destroy()
    call thermal%destroy()
    call river%destroy()
  end subroutine check_accumulation_cadence

  !> \brief Check mean flow/lake publication, buffer swaps, and stable exchange targets over two substeps.
  subroutine check_flow_completion()
    type(grid_t), target :: grid
    type(grid_t), pointer :: grid_ptr
    type(river_t), target :: river
    type(river_flow_t) :: processor
    type(river_flow_state_t), target :: state
    type(river_muskingum_state_t) :: muskingum
    real(dp), pointer :: discharge_target(:), lake_target(:), input_buffer(:), routed_buffer(:)
    integer(i8) :: n_levels

    call grid%init(3_i4,1_i4,cellsize=2000.0_dp)
    grid_ptr => grid
    call river%from_fdir([d8_E,d8_E,sink],grid_ptr)
    call river%calc_order(root=.false.,sort=.true.)
    n_levels = river%order%n_levels
    call state%allocate(river%n_nodes,1_i4)
    call processor%init(river,3_i8,3_i8,.false.,.false.,lake_nodes=[2_i8])
    call muskingum%allocate(river%n_nodes)
    muskingum%nu1 = 1.0_dp
    muskingum%nu2 = 1.0_dp
    muskingum%iterations = 2_i4
    discharge_target => state%discharge
    lake_target => state%lake_inflow
    input_buffer => state%current_input
    routed_buffer => state%current_routed

    processor%runoff = [2.0_dp,3.0_dp,7.0_dp]
    call processor%route_nodes(1_i8,n_levels,state,muskingum,lake_outflow=[4.0_dp])
    call assert_true(all(state%current_input==[2.0_dp,4.0_dp,11.0_dp]),'lake release and raw node flow')
    call processor%advance(state)
    call assert_true(associated(input_buffer,state%previous_input),'flow input history swaps without a copy')
    call assert_true(associated(routed_buffer,state%previous_routed),'flow routed history swaps without a copy')
    call assert_true(all(state%discharge==0.0_dp) .and. all(state%lake_inflow==0.0_dp),'substeps leave published flow means unchanged')
    processor%runoff = [4.0_dp,6.0_dp,9.0_dp]
    call processor%route_nodes(1_i8,n_levels,state,muskingum,lake_outflow=[8.0_dp])
    call processor%advance(state)
    call processor%complete(state,muskingum)
    call assert_true(all(state%discharge==[3.0_dp,6.0_dp,14.0_dp]),'completed node discharge averages all substeps')
    call assert_true(all(state%lake_inflow==7.5_dp),'canonical lake inflow averages all substeps')
    call assert_true(associated(discharge_target,state%discharge) .and. associated(lake_target,state%lake_inflow),'completion preserves published targets')
    call state%reset()
    call assert_true(all(state%discharge==0.0_dp) .and. all(state%lake_inflow==0.0_dp),'state reset clears completed flow means')
    call assert_true(associated(discharge_target,state%discharge) .and. associated(lake_target,state%lake_inflow),'reset preserves published targets')
    call processor%destroy()
    call assert_true(allocated(state%lake_inflow) .and. allocated(state%current_input),'processor destruction preserves owned flow state')
    nullify(discharge_target,lake_target,input_buffer,routed_buffer)
    call state%destroy()
    call muskingum%destroy()
    call river%destroy()
  end subroutine check_flow_completion

  !> \brief Check router-owned land-cell mapping, flow use, and rebuilding with different runoff support.
  subroutine check_shared_node_map()
    type(grid_t), target :: full, land
    type(grid_t), pointer :: full_ptr, land_ptr
    type(river_t), target :: river
    type(river_t), pointer :: river_ptr
    type(river_router_t) :: router
    type(river_flow_state_t) :: flow
    type(river_thermal_state_t) :: thermal
    real(dp) :: expected(3)
    logical :: completed

    call full%init(3_i4,1_i4,cellsize=2000.0_dp)
    call land%init(3_i4,1_i4,cellsize=2000.0_dp,mask=reshape([.true.,.false.,.true.],[3,1]))
    full_ptr => full
    land_ptr => land
    call river%from_fdir([d8_E,d8_E,sink],full_ptr)
    ! Represent one SCC node per full-grid cell, with the middle cell outside land support.
    river%scc = .true.
    allocate(river%node_cell(3),source=[1_i8,2_i8,3_i8])
    allocate(river%area_fraction(3),source=1.0_dp)
    river_ptr => river
    call flow%allocate(river%n_nodes)
    call router%init(river_ptr,flow,thermal,[1.0_dp,1.0_dp,1.0_dp],land_ptr,runoff_grid=land_ptr,omp_level_thresh=0_i8)
    call assert_true(all(router%node_runoff_cell==[1_i8,0_i8,2_i8]),'router owns full-to-land node map')
    call router%update([1.0_dp,2.0_dp],flow,thermal,completed)
    expected = [1.0_dp,0.0_dp,2.0_dp]*(4000000.0_dp*0.001_dp/3600.0_dp)
    call assert_true(all(abs(router%flow%runoff-expected)<1.e-12_dp),'flow uses shared map and omits masked runoff cells')
    call router%flow%destroy()
    call assert_true(allocated(router%node_runoff_cell),'processor destruction preserves shared map')
    call router%init(river_ptr,flow,thermal,[1.0_dp,1.0_dp,1.0_dp],full_ptr,omp_level_thresh=0_i8)
    call assert_true(.not.allocated(router%node_runoff_cell),'same-grid reinitialization removes obsolete node map')
    call router%destroy()
    call flow%destroy()
    call river%destroy()
  end subroutine check_shared_node_map

  !> \brief Read updated external sealed-fraction slices during temperature-input preparation.
  subroutine check_sealed_fraction()
    type(domain_t), target :: domain
    real(dp), allocatable, target :: sealed_first(:), sealed_second(:), baseflow(:), sealed_runoff(:)
    real(dp), allocatable :: expected(:), air(:)
    integer(i4) :: step
    integer(i8) :: n
    call domain%create(main_file="test_nml/mrm_rivertemp_meteo_minimal.nml",cwd=".")
    call domain%configure(main_file="test_nml/mrm_rivertemp_meteo_minimal.nml", &
      para_file="test_nml/mhm_parameter_v6_meteo.nml")
    call domain%exchange%f_sealed%provide("test")
    call domain%exchange%baseflow%provide("test")
    call domain%exchange%runoff_sealed%provide("test")
    call domain%connect()
    n = domain%exchange%level1_land%ncells
    allocate(sealed_first(n),source=0.1_dp)
    allocate(sealed_second(n),source=0.7_dp)
    allocate(baseflow(n),source=0.01_dp)
    allocate(sealed_runoff(n),source=0.02_dp)
    allocate(expected(n),air(n))
    call domain%exchange%f_sealed%prepare_data("test",24_i4)
    call domain%exchange%baseflow%prepare_data("test",1_i4)
    call domain%exchange%runoff_sealed%prepare_data("test",1_i4)
    domain%exchange%f_sealed%data => sealed_first
    domain%exchange%baseflow%data => baseflow
    domain%exchange%runoff_sealed%data => sealed_runoff
    call domain%initialize()
    do step = 1_i4, 2_i4
      if (step == 2_i4) domain%exchange%f_sealed%data => sealed_second
      call domain%update()
      air = max(T0_dp,domain%exchange%temp%data+T0_dp)
      expected = domain%exchange%runoff_total%data*air + &
        baseflow*(1.0_dp-domain%exchange%f_sealed%data)* &
        (max(T0_dp+5.0_dp,domain%exchange%tann%data+T0_dp)-air) + &
        sealed_runoff*domain%exchange%f_sealed%data*(max(T0_dp,domain%exchange%temp%data+T0_dp-1.5_dp)-air)
      call assert_true(all(abs(expected-domain%mrm%temperature_input%depth)<1.e-10_dp),'temporal sealed-fraction input')
    end do
    call domain%destroy()
  end subroutine check_sealed_fraction

  !> \brief Fail with a descriptive invariant instead of accepting a silent numerical difference.
  subroutine assert_true(condition,message)
    logical, intent(in) :: condition !< Invariant under test.
    character(*), intent(in) :: message !< Failure diagnostic.
    if (.not.condition) then
      print *, message
      error stop 1
    end if
  end subroutine assert_true
end program test_river_routing_optimization
