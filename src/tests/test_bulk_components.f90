!> \file test_bulk_components.f90
!> \brief Exercise meteo PET and output preparation with a substantial parallel workload.
!> \date Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
program test_bulk_components
  use mo_kind, only: dp, i4, i8
  use mo_exchange_type, only: exchange_t, var_dp
  use mo_grid, only: grid_t, spherical
  use mo_datetime, only: datetime
  use mo_grid_io, only: daily
  use mo_meteo_container, only: meteo_t
  use mo_mhm_container, only: mhm_t
  use mo_pet, only: pet_hargreaves, pet_priestly, pet_penman
  use mo_meteo_constants, only: HarSamConst
  implicit none
  integer, parameter :: n = 65536
  type(exchange_t), target :: exchange
  type(grid_t), target :: grid
  type(meteo_t), target :: meteo
  type(mhm_t), target :: mhm
  real(dp), allocatable, target :: temp(:), tmin(:), tmax(:), netrad(:), eabs(:), wind(:), factor(:)
  real(dp), allocatable, target :: moisture(:,:), saturation(:,:), infiltration(:,:), aet(:,:), sealed(:)
  real(dp) :: expected
  integer(i8) :: start, finish, rate
  integer(i4) :: status, selection, j, k
  integer(i4), parameter :: cases(5) = [-2_i4,-1_i4,1_i4,2_i4,3_i4]
  allocate(temp(n),tmin(n),tmax(n),netrad(n),eabs(n),wind(n),factor(n))
  temp = 10.0_dp
  tmin = 5.0_dp
  tmax = 15.0_dp
  netrad = 100.0_dp
  eabs = 1000.0_dp
  wind = 2.0_dp
  factor = 1.0_dp
  do selection = 1, size(cases)
    call setup_exchange()
    status = exchange%config%processes%set(pet=cases(selection))
    if (status /= 0) error stop 'process configuration failed'
    call provide(exchange%raw_pre,factor,daily)
    call provide(exchange%raw_temp,temp,daily)
    call provide(exchange%raw_tmin,tmin,daily)
    call provide(exchange%raw_tmax,tmax,daily)
    call provide(exchange%raw_netrad,netrad,daily)
    call provide(exchange%raw_eabs,eabs,daily)
    call provide(exchange%raw_wind,wind,daily)
    call provide(exchange%raw_pet,factor,daily)
    call provide(exchange%pet_fac_aspect,factor,0_i4)
    call provide(exchange%pet_fac_lai,factor,0_i4)
    call provide(exchange%pet_coeff_hs,factor,0_i4)
    call provide(exchange%pet_coeff_pt,factor,0_i4)
    call provide(exchange%resist_aero,factor,0_i4)
    call provide(exchange%resist_surf,factor,0_i4)
    status = meteo%config%set()
    if (status /= 0) error stop 'meteo configuration failed'
    meteo%exchange => exchange
    call meteo%configure()
    call meteo%connect()
    call meteo%initialize()
    call system_clock(start,rate)
    do j = 1, 20
      call meteo%update()
    end do
    call system_clock(finish)
    write(*,'(a,i0,a,f10.4)') 'bulk PET case ',cases(selection),' seconds: ',real(finish-start,dp)/real(rate,dp)
    do k = 1, n
      select case(cases(selection))
      case(-2,-1)
        expected = 1.0_dp
      case(1)
        expected = pet_hargreaves(1.0_dp,HarSamConst,10.0_dp,15.0_dp,5.0_dp, &
          meteo%scratch%latitude(k),exchange%time_step_start%doy())
      case(2)
        expected = pet_priestly(1.0_dp,100.0_dp,10.0_dp)
      case(3)
        expected = pet_penman(100.0_dp,10.0_dp,1.0_dp,0.5_dp,1.0_dp,1.0_dp,1.0_dp)
      end select
      if (abs(exchange%pet%data(k)-expected) > 1.e-12_dp * max(1.0_dp,abs(expected))) error stop 'bulk PET changed values'
    end do
    call exchange%destroy()
    call meteo%destroy()
    call grid%destroy()
  end do
  call setup_exchange()
  status = exchange%config%processes%set(pet=0_i4)
  if (status /= 0) error stop 'process reset failed'
  exchange%n_layers = 2
  allocate(moisture(n,2),saturation(n,2),infiltration(n,2),aet(n,2),sealed(n))
  moisture = 30.0_dp
  saturation = 100.0_dp
  infiltration = 2.0_dp
  aet = 1.0_dp
  sealed = 0.2_dp
  call exchange%soil_moisture%provide('test')
  call exchange%soil_moisture%prepare_data('test',24_i4)
  exchange%soil_moisture%data => moisture
  call exchange%sm_saturation%provide('test')
  call exchange%sm_saturation%prepare_data('test',24_i4)
  exchange%sm_saturation%data => saturation
  call exchange%infiltration%provide('test')
  call exchange%infiltration%prepare_data('test',24_i4)
  exchange%infiltration%data => infiltration
  call exchange%aet_soil%provide('test')
  call exchange%aet_soil%prepare_data('test',24_i4)
  exchange%aet_soil%data => aet
  call provide(exchange%f_sealed,sealed,0_i4)
  mhm%exchange => exchange
  mhm%soil%horizon_bounds = [0.0_dp,200.0_dp,1000.0_dp]
  mhm%soil%variable_bottom = .false.
  exchange%model_horizon_bounds => mhm%soil%horizon_bounds
  exchange%model_horizon_bottom_variable = .false.
  mhm%io%output_active = .true.
  mhm%io%calc_f_not_sealed = .true.
  mhm%io%output_path = 'bulk_output.nc'
  status = mhm%config%set()
  if (status /= 0) error stop 'mHM configuration failed'
  status = mhm%output_config%set(output_frequency=0_i4, output_deflate_level=0_i4, &
    out_swc=.true.,out_sm=.true.,out_sm_all=.true.,out_soil_infil=.true.,out_aet_layer=.true.)
  if (status /= 0) error stop 'output configuration failed'
  call mhm%initialize()
  call system_clock(start,rate)
  do j = 1, 20
    call mhm%update()
    if (any(mhm%output_work%pervious /= 1.0_dp-sealed)) error stop 'stale pervious fraction'
    sealed = sealed + 0.001_dp
  end do
  call system_clock(finish)
  write(*,'(a,f10.4)') 'bulk output preparation seconds: ',real(finish-start,dp)/real(rate,dp)
  call mhm%finalize()
  call exchange%destroy()
  call mhm%destroy()
  call grid%destroy()
contains
  !> \brief Provide real container inputs on stable, application-owned storage.
  subroutine provide(variable, values, stepping)
    type(var_dp), intent(inout) :: variable !< Exchange field.
    real(dp), intent(in), target :: values(:) !< Stable input storage.
    integer(i4), intent(in) :: stepping !< Native timing.
    integer(i4) :: native
    native = stepping
    if (.not.variable%static .and. stepping == 0_i4) native = 24_i4
    call variable%provide('test')
    call variable%prepare_data('test',native)
    variable%data => values
  end subroutine provide

  !> \brief Set up a large workload for the component parallel loops.
  subroutine setup_exchange()
    status = exchange%config%project%set()
    if (status /= 0) error stop 'project configuration failed'
    call exchange%create(cwd='.')
    exchange%step_hours = 24
    exchange%time_step_start = datetime('2001-01-02 00:00:00')
    exchange%time = datetime('2001-01-03 00:00:00')
    exchange%start_time = exchange%time_step_start
    exchange%end_time = datetime('2001-02-01 00:00:00')
    call grid%init(256_i4,256_i4,yllcorner=40.0_dp,cellsize=0.01_dp,coordsys=spherical)
    exchange%level0_land => grid
    exchange%level0 => grid
    exchange%level1_land => grid
    exchange%level2 => grid
  end subroutine setup_exchange
end program test_bulk_components
