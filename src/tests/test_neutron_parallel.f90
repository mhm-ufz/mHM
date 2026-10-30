!> \file test_neutron_parallel.f90
!> \brief Shared neutron geometry remains usable across initialization and OpenMP team changes.
!> \date Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
program test_neutron_parallel
  use mo_kind, only: dp, i4
  use mo_exchange_type, only: exchange_t, var_dp
  use mo_grid, only: grid_t, cartesian
  use mo_datetime, only: datetime
  use mo_mhm_container, only: mhm_t
  use mo_neutrons, only: DesiletsN0, COSMIC
#ifdef _OPENMP
  use omp_lib, only: omp_set_num_threads
#endif
  implicit none
  integer(i4), parameter :: n = 65536, layers = 4
  type(exchange_t), target :: exchange
  type(grid_t), target :: grid
  type(mhm_t), target :: mhm
  real(dp), target :: bounds(layers+1) = [0.0_dp,103.5_dp,307.2_dp,1005.5_dp,1500.0_dp]
  real(dp), allocatable, target :: moisture(:,:), density(:,:), lattice(:,:), l3(:,:), n0(:), interception(:), snow(:)
  real(dp), allocatable :: reference(:)
  integer(i4) :: neutron_case, status, k, repetition, step

  allocate(moisture(n,layers), density(n,layers), lattice(n,layers), l3(n,layers), n0(n), interception(n), snow(n), reference(n))
  density = 1.4_dp
  lattice = 2.0_dp
  l3 = 90.0_dp
  n0 = 1500.0_dp
  do neutron_case = 1, 2
    status = exchange%config%project%set()
    if (status /= 0) error stop 'project configuration failed'
    call exchange%create(cwd='.')
    status = exchange%config%processes%set(neutrons=neutron_case)
    if (status /= 0) error stop 'neutron configuration failed'
    exchange%n_layers = layers
    exchange%step_hours = 24
    exchange%time_step_start = datetime('2001-01-02 00:00:00')
    exchange%time = datetime('2001-01-03 00:00:00')
    exchange%start_time = exchange%time_step_start
    exchange%end_time = exchange%time
    call grid%init(256_i4,256_i4,cellsize=1000.0_dp,coordsys=cartesian)
    exchange%level1_land => grid
    exchange%model_horizon_bounds => bounds
    exchange%model_horizon_bottom_variable = .false.
    call exchange%soil_moisture%provide('test')
    call exchange%soil_moisture%prepare_data('test',24_i4)
    exchange%soil_moisture%data => moisture
    call exchange%bulk_density%provide('test')
    exchange%bulk_density%static = .true.
    call exchange%bulk_density%prepare_data('test',0_i4)
    exchange%bulk_density%data => density
    call exchange%lattice_water%provide('test')
    exchange%lattice_water%static = .true.
    call exchange%lattice_water%prepare_data('test',0_i4)
    exchange%lattice_water%data => lattice
    call exchange%cosmic_l3%provide('test')
    exchange%cosmic_l3%static = .true.
    call exchange%cosmic_l3%prepare_data('test',0_i4)
    exchange%cosmic_l3%data => l3
    call provide(exchange%desilets_n0, n0, 0_i4)
    call provide(exchange%interception, interception, 24_i4)
    call provide(exchange%snowpack, snow, 24_i4)
    mhm%exchange => exchange
    status = mhm%config%set()
    if (status /= 0) error stop 'mHM configuration failed'
    status = mhm%output_config%init()
    if (status /= 0) error stop 'mHM output defaults failed'
    call mhm%configure()
    call mhm%prepare_restart()
    call mhm%connect()
    do repetition = 1, 2
#ifdef _OPENMP
      call omp_set_num_threads(1)
#endif
      call mhm%initialize()
      do step = 1, 2
        do k = 1, n
          moisture(k,1) = 15.0_dp + real(mod(k,17),dp) * real(repetition*step,dp)
          moisture(k,2) = 40.0_dp + real(mod(k,23),dp) * real(repetition*step,dp)
          moisture(k,3) = 180.0_dp + real(mod(k,31),dp) * real(repetition*step,dp)
          moisture(k,4) = 800.0_dp ! Excluded bottom layer, deliberately distinct.
          interception(k) = real(mod(k,3),dp) * real(step,dp)
          snow(k) = real(mod(k,5),dp) * real(repetition,dp)
          if (neutron_case == 1) then
            call DesiletsN0(moisture(k,:layers-1),bounds(2:layers),density(k,:layers-1),lattice(k,:layers-1),n0(k),reference(k))
          else
            call COSMIC(moisture(k,:layers-1),bounds(2:layers),mhm%neutrons%integral_afast,interception(k),snow(k), &
              n0(k),density(k,:layers-1),lattice(k,:layers-1),l3(k,:layers-1),reference(k))
          end if
        end do
#ifdef _OPENMP
        if (step == 1) call omp_set_num_threads(4)
        if (step == 2) call omp_set_num_threads(2)
#endif
        call mhm%update()
        if (any(abs(exchange%neutrons%data-reference) > 1.e-12_dp*max(1.0_dp,abs(reference)))) &
          error stop 'shared geometry or changed OpenMP team altered neutron results'
      end do
    end do
    call mhm%finalize()
    call exchange%destroy()
    call mhm%destroy()
    call grid%destroy()
  end do
contains
  !> \brief Publish stable, externally owned input storage.
  subroutine provide(variable, values, stepping)
    type(var_dp), intent(inout) :: variable !< Exchange field.
    real(dp), intent(in), target :: values(:) !< Stable application-owned values.
    integer(i4), intent(in) :: stepping !< Native timing of the published field.
    call variable%provide('test')
    variable%static = stepping == 0_i4
    call variable%prepare_data('test',stepping)
    variable%data => values
  end subroutine provide
end program test_neutron_parallel
