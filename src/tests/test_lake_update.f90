program test_lake_update
  use mo_kind, only: dp, i4, i8
  use mo_datetime, only: datetime
  use mo_exchange_type, only: exchange_t
  use mo_mlm_container, only: mlm_t
  use mo_lake_bathymetry, only: generate_lake_bathymetry
  implicit none
  integer, parameter :: n = 8
  type(exchange_t), target :: exchange
  type(mlm_t), target :: mlm
  real(dp), target :: inflow(n), pre(n), air(n), sw(n), lw(n)
  real(dp) :: before(n), residual(n), evaporation_volume(n), footprint(n)
  integer :: status, i, phase
  status = exchange%config%project%set(max_lakes=n)
  call check(status == 0, 'project configuration')
  status = exchange%config%processes%set(lake=1_i4,routing=2_i4)
  call check(status == 0, 'process configuration')
  call exchange%create(cwd='.')
  exchange%step_hours = 1_i4
  exchange%time_step_start = datetime('2001-06-01 12:00:00')
  exchange%time = datetime('2001-06-01 13:00:00')
  mlm%exchange => exchange
  call mlm%set_dims()
  status = mlm%config%set(evaporation_method=2_i4, varying_albedo=.false., &
    varying_lambda=.false., varying_psychro=.false.)
  call check(status == 0, 'lake configuration')
  mlm%definition%ids = [(int(i,i8),i=1,n)]
  mlm%definition%spill_reference_level = spread(100.0_dp,1,n)
  mlm%definition%spill_reference_surface_area = spread(1.0e4_dp,1,n)
  mlm%definition%spill_reference_volume = spread(1.0e4_dp,1,n)
  mlm%definition%footprint_area = spread(1.0e4_dp,1,n)
  mlm%definition%latitude = spread(48.0_dp,1,n)
  mlm%definition%longitude = spread(11.0_dp,1,n)
  mlm%definition%ks = spread(0.0_dp,1,n)
  mlm%definition%kbf = spread(100.0_dp,1,n)
  mlm%definition%spillway_length = spread(0.0_dp,1,n)
  mlm%definition%endorheic = spread(.false.,1,n)
  mlm%runtime%b = 1.0_dp
  mlm%runtime%cd = 1.8_dp
  mlm%runtime%exponent = 2.0_dp
  allocate(mlm%state%hav(n),mlm%state%level(n),mlm%state%surface_area(n),mlm%state%volume(n),mlm%state%gw_storage(n))
  allocate(mlm%state%level_previous(n),mlm%state%surface_area_previous(n), &
    mlm%state%volume_previous(n),mlm%state%gw_storage_previous(n))
  allocate(mlm%state%outflow(n),mlm%state%evaporation(n),mlm%state%seepage(n),mlm%state%spill(n),mlm%state%baseflow(n))
  do i = 1,n
    call generate_lake_bathymetry(1_i4,100.0_dp,1.0e4_dp,1.0e4_dp, &
      mlm%state%hav(i)%level,mlm%state%hav(i)%surface_area,mlm%state%hav(i)%volume)
  end do
  inflow = 0.0_dp
  pre = 0.0_dp
  air = 20.0_dp
  sw = 1.0e8_dp
  lw = 350.0_dp
  call exchange%lake_inflow%provide('test')
  call exchange%lake_inflow%prepare_data('test',1_i4)
  exchange%lake_inflow%data => inflow
  call exchange%lake_pre%provide('test')
  call exchange%lake_pre%prepare_data('test',1_i4)
  exchange%lake_pre%data => pre
  call exchange%lake_temp%provide('test')
  call exchange%lake_temp%prepare_data('test',1_i4)
  exchange%lake_temp%data => air
  call exchange%lake_ssrd%provide('test')
  call exchange%lake_ssrd%prepare_data('test',1_i4)
  exchange%lake_ssrd%data => sw
  call exchange%lake_strd%provide('test')
  call exchange%lake_strd%prepare_data('test',1_i4)
  exchange%lake_strd%data => lw
  do phase = 1,3
    ! Constant-area bathymetry: bottom=99 m; h=101 m gives exactly 20,000 m3.
    mlm%state%level = 101.0_dp
    mlm%state%surface_area = 1.0e4_dp
    mlm%state%volume = 2.0e4_dp
    mlm%state%gw_storage = 0.0_dp
    mlm%state%outflow = 0.0_dp
    mlm%definition%ks = 0.0_dp
    mlm%definition%spillway_length = 0.0_dp
    mlm%definition%endorheic = .false.
    mlm%evaporation_active = phase == 2
    if (phase == 1) then
      mlm%definition%spillway_length(1:4) = 1.0e3_dp
      mlm%definition%ks(5:8) = 1.0e6_dp
    else if (phase == 2) then
      ! All three demanded losses exhaust storage. Include closed lakes in the same team.
      mlm%definition%spillway_length(1:4) = 1.0e3_dp
      mlm%definition%ks(1:4) = 1.0e6_dp
      mlm%definition%endorheic(5:8) = .true.
    else
      mlm%state%level = 100.0_dp
      mlm%state%volume = 1.0e4_dp
      footprint = [5000.0_dp,15000.0_dp,5000.0_dp,15000.0_dp,5000.0_dp,15000.0_dp,5000.0_dp,15000.0_dp]
      mlm%definition%footprint_area = footprint
      pre = 10.0_dp
    end if
    mlm%state%level_previous = mlm%state%level
    mlm%state%surface_area_previous = mlm%state%surface_area
    mlm%state%volume_previous = mlm%state%volume
    mlm%state%gw_storage_previous = mlm%state%gw_storage
    before = mlm%state%volume
    call mlm%update()
    if (phase < 3) then
      evaporation_volume = mlm%state%evaporation/1000.0_dp * 0.5_dp*(1.0e4_dp+mlm%state%surface_area)
      residual = before - mlm%state%volume - mlm%state%gw_storage - mlm%state%outflow*3600.0_dp - evaporation_volume
      call check(maxval(abs(residual)) < 1.0e-6_dp, 'drying does not conserve total lake + groundwater + routed + evaporated water')
      call check(maxval(abs(mlm%state%volume)) < 1.0e-6_dp, 'exhausted lake retained volume')
      if (phase == 1) then
        call check(maxval(abs(mlm%state%spill(1:4)*3600.0_dp-before(1:4))) < 1.0e-6_dp, 'spill was not capped')
        call check(maxval(abs(mlm%state%gw_storage(5:8)-before(5:8))) < 1.0e-6_dp, 'seepage was not capped before groundwater transfer')
      else
        call check(maxval(abs(mlm%state%spill)) < 1.0e-6_dp, 'spill consumed unavailable water before evaporation')
        call check(maxval(abs(mlm%state%seepage)) < 1.0e-6_dp, 'seepage consumed unavailable water before evaporation')
        call check(maxval(abs(evaporation_volume-before)) < 1.0e-6_dp, 'evaporation was not capped')
      end if
    else
      ! P alone must add P * footprint, for both exposed and inundated bank rings.
      call check(maxval(abs(mlm%state%volume-before-pre/1000.0_dp*footprint)) < 1.0e-6_dp, 'bank-area signs violate the footprint precipitation budget')
    end if
  end do
  call exchange%destroy()
  call mlm%destroy()
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ', message
      error stop 1
    end if
  end subroutine check
end program test_lake_update
