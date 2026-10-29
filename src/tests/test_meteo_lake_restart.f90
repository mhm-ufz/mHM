program test_meteo_lake_restart
  use mo_kind, only: sp, dp, i4, i8
  use mo_datetime, only: datetime
  use mo_exchange_type, only: exchange_t
  use mo_grid, only: grid_t, bottom_up
  use mo_points, only: points_t, cartesian
  use mo_grid_io, only: no_time, daily
  use mo_lake_forcing, only: lake_support_build
  use mo_meteo_lake_restart, only: lake_forcing_restart_t, lake_grids_match
  use mo_meteo_container, only: meteo_t
  use mo_netcdf, only: NcDataset, NcVariable, NcDimension
  implicit none
  type(lake_forcing_restart_t), target :: cache, restored
  type(exchange_t), target :: exchange
  type(meteo_t), target :: meteo
  type(grid_t), target :: level2, geometry
  type(points_t), target :: points
  integer(i8), target :: ids(3), map(6)
  real(dp), target :: raw(6)
  real(dp), pointer :: published(:)
  real(dp) :: weight, corrupt_weights(6)
  type(NcDataset) :: nc
  type(NcVariable) :: var
  type(NcDimension) :: dimensions(4)
  real(dp) :: source_weights(2,2,12,24)
  integer :: status, unit
  logical :: fractions, rebuild
  character(64) :: mode
  call get_command_argument(1, mode)
  fractions = mode == 'fractions' .or. mode == 'rebuild' .or. mode == 'default_mode' .or. &
    mode == 'other_domain' .or. mode == 'missing_fraction' .or. mode == 'fraction_conflict'
  rebuild = mode == 'rebuild' .or. mode == 'rebuild_weights' .or. mode == 'rebuild_weighted' .or. mode == 'weight_origin'
  ids = [101_i8,205_i8,309_i8]
  map = [101_i8,101_i8,205_i8,309_i8,309_i8,309_i8]
  call level2%init(3_i4,2_i4, cellsize=2.0_dp, coordsys=cartesian, y_direction=bottom_up)
  geometry = level2
  cache%grid = level2
  cache%ids = ids
  call lake_support_build(level2, geometry, map, ids, cache%support)
  cache%weighted = .true.
  cache%disaggregate_pre = .true.
  cache%disaggregate_pet = .true.
  allocate(cache%pre_weights(3,12,24), cache%pet_weights(3,12,24))
  weight = 1.0_dp/24.0_dp
  if (mode == 'floatweights') weight = real(real(weight,sp),dp)
  cache%pre_weights = weight
  cache%pet_weights = weight
  if (fractions) then
    cache%weighted = .false.
    cache%frac_pre = 0.25_dp
    cache%frac_pet = 0.1_dp
    deallocate(cache%pre_weights,cache%pet_weights)
  end if
  call cache%write('support.nc', '1990-01-01 00:00:00')
  select case(trim(mode))
  case('empty','negative','nan','bounds','duplicate','normalization','temporal','offsets')
    ! Corrupt a valid file so each failure exercises the production reader, not its writer.
    nc = NcDataset('support.nc','a')
    corrupt_weights = [0.5_dp,0.5_dp,1.0_dp,1.0_dp/3.0_dp,1.0_dp/3.0_dp,1.0_dp/3.0_dp]
    select case(trim(mode))
    case('empty')
      var = nc%getVariable('support_count')
      call var%setData([0_i8,1_i8,3_i8])
    case('negative')
      var = nc%getVariable('support_weight')
      corrupt_weights(1:2) = [-1.0_dp,2.0_dp]
      call var%setData(corrupt_weights)
    case('nan')
      block
        use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
        var = nc%getVariable('support_weight')
        corrupt_weights(1) = ieee_value(0.0_dp,ieee_quiet_nan)
        call var%setData(corrupt_weights)
      end block
    case('bounds')
      var = nc%getVariable('support_l2_index')
      call var%setData([0_i8,2_i8,3_i8,4_i8,5_i8,6_i8])
    case('duplicate')
      var = nc%getVariable('lake_id')
      call var%setData([101_i8,101_i8,309_i8])
    case('normalization')
      var = nc%getVariable('support_weight')
      corrupt_weights(1:2) = 1.0_dp
      call var%setData(corrupt_weights)
    case('temporal')
      var = nc%getVariable('pre_weights')
      call var%setData(cache%pre_weights*24.0_dp)
    case('offsets')
      var = nc%getVariable('support_start')
      call var%setData([1_i8,1_i8,4_i8])
    end select
    call nc%close()
  end select
  call restored%read('support.nc')
  call check(lake_grids_match(restored%grid,level2), 'grid definition did not round trip')
  nc = NcDataset('support.nc', 'r')
  call check(.not.nc%hasVariable('L2_area'), 'support restart redundantly writes cell areas')
  call nc%close()
  status = exchange%config%project%set(n_domains=2_i4)
  call check(status == 0, 'project configuration')
  status = exchange%config%processes%set(lake=-2_i4, routing=2_i4)
  call check(status == 0, 'process configuration')
  call exchange%create(cwd='.')
  exchange%cwd = exchange%root
  exchange%step_hours = 1_i4
  exchange%time_step_start = datetime('2007-06-01 12:00:00')
  exchange%time = datetime('2007-06-01 13:00:00')
  exchange%level2 => level2
  call points%init([1.0_dp,3.0_dp,5.0_dp], [1.0_dp,1.0_dp,3.0_dp], coordsys=cartesian)
  exchange%lake_points => points
  ids = [205_i8,101_i8,309_i8]
  if (mode == 'ids') ids(1) = 999_i8
  call exchange%lake_ids%provide('test')
  call exchange%lake_ids%prepare_data('test', no_time)
  exchange%lake_ids%data => ids
  raw = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
  call exchange%raw_pre%provide('Input')
  call exchange%raw_pre%prepare_data('Input', daily)
  exchange%raw_pre%data => raw
  call exchange%raw_pet%provide('Input')
  call exchange%raw_pet%prepare_data('Input', daily)
  exchange%raw_pet%data => raw
  meteo%exchange => exchange
  call meteo%set_dims()
  status = meteo%config%set(read_restart=[.true.], restart_input_path=[character(256)::'support.nc'], &
    write_restart=[.true.], restart_output_path=[character(256)::'support_reordered.nc'], read_meteo_weights=[.not.fractions])
  call check(status == 0, 'meteo configuration')
  if (mode == 'conflict') meteo%config%read_meteo_weights(1) = .false.
  if (fractions .or. mode == 'omitted_conflict' .or. mode == 'explicit_false' .or. mode == 'file_weighted') then
    open(newunit=unit,file='meteo_support.nml',status='replace')
    write(unit,'(a)') '&config_meteo', &
      'read_restart(1)=.true., restart_input_path(1)="support.nc",', &
      'write_restart(1)=.true., restart_output_path(1)="support_reordered.nc",'
    if (mode == 'explicit_false') write(unit,'(a)') 'READ_METEO_WEIGHTS(1)=.false., ! explicit conflict'
    if (mode == 'other_domain') write(unit,'(a)') 'read_meteo_weights(2)=.true., ! domain 1 uses its false default'
    if (mode == 'file_weighted') write(unit,'(a)') 'read_meteo_weights(1)=.true.,'
    if (fractions) then
      if (mode /= 'missing_fraction') write(unit,'(a)') 'frac_night_pre(:,1)=12*0.25,'
      if (mode == 'fraction_conflict') write(unit,'(a)') 'frac_night_pre(:,1)=12*0.5,'
      write(unit,'(a)') 'frac_night_pet(:,1)=12*0.1,'
    end if
    write(unit,'(a)') '/'
    close(unit)
    call meteo%configure(file='meteo_support.nml')
  else
    call meteo%configure()
  end if
  call meteo%prepare_restart()
  ! Input-owned geometry must survive validation. Test missing and changed geometry as well.
  if (mode == 'geometry' .or. mode == 'mapping' .or. rebuild) then
    exchange%level0_lake => geometry
    call exchange%lake_map%provide('Input')
    call exchange%lake_map%prepare_data('Input', no_time)
    exchange%lake_map%data => map
    if (mode == 'mapping') map(1) = 205_i8
  end if
  if (rebuild) then
    call level2%destroy()
    call level2%init(2_i4,2_i4,cellsize=4.0_dp,coordsys=cartesian,y_direction=bottom_up)
    exchange%raw_pre%data => raw(:4)
    exchange%raw_pet%data => raw(:4)
    if (mode == 'rebuild_weighted' .or. mode == 'weight_origin') then
      nc = NcDataset('sourceweights.nc','w')
      if (mode == 'weight_origin') level2%xllcorner = 1.0_dp
      call level2%to_netcdf(nc)
      if (mode == 'weight_origin') level2%xllcorner = 0.0_dp
      dimensions(1) = nc%getDimension('x')
      dimensions(2) = nc%getDimension('y')
      dimensions(3) = nc%setDimension('month',12)
      dimensions(4) = nc%setDimension('hour',24)
      source_weights = 1.0_dp/24.0_dp
      var = nc%setVariable('pre_weights','f64',dimensions)
      call var%setAttribute('missing_value',-9999.0_dp)
      call var%setData(source_weights)
      var = nc%setVariable('pet_weights','f64',dimensions)
      call var%setAttribute('missing_value',-9999.0_dp)
      call var%setData(source_weights)
      call nc%close()
      meteo%config%pre_weights_path(1) = 'sourceweights.nc'
      meteo%config%pet_weights_path(1) = 'sourceweights.nc'
    end if
  end if
  if (mode == 'origin') level2%xllcorner = 1.0_dp
  if (mode == 'resolution') level2%cellsize = 3.0_dp
  if (mode == 'orientation') level2%y_direction = 0_i4
  if (mode == 'mask') then
    call level2%destroy()
    call level2%init(3_i4,2_i4, cellsize=2.0_dp, coordsys=cartesian, y_direction=bottom_up, &
      mask=reshape([.true.,.true.,.false.,.true.,.true.,.true.],[3,2]))
  end if
  call meteo%connect()
  if (mode == 'geometry') call check(associated(exchange%level0_lake,geometry), 'replaced Input-owned geometry')
  published => exchange%lake_pre%data
  call meteo%initialize()
  call meteo%update()
  if (mode == 'rebuild_weighted') then
    call check(maxval(abs(exchange%lake_pre%data-[2.0_dp,1.0_dp,4.0_dp/3.0_dp]/24.0_dp)) < 1.0e-13_dp, &
      'changed L2 did not rebuild weighted support from compatible source files')
  else if (mode == 'rebuild') then
    call check(maxval(abs(exchange%lake_pre%data-[2.0_dp,1.0_dp,4.0_dp/3.0_dp]*1.5_dp/24.0_dp)) < 1.0e-13_dp, &
      'changed L2 did not rebuild from supplied geometry')
  else if (fractions) then
    call check(maxval(abs(exchange%lake_pre%data-[3.0_dp,1.5_dp,5.0_dp]*1.5_dp/24.0_dp)) < 1.0e-13_dp, &
      'configured day/night fractions did not match saved support')
  else
    call check(maxval(abs(exchange%lake_pre%data-[3.0_dp,1.5_dp,5.0_dp]*weight)) < 1.0e-13_dp, &
      'cache did not reorder stable IDs or reuse across dates')
  end if
  call meteo%initialize()
  call check(associated(published,exchange%lake_pre%data), 'reinitialization replaced forcing target')
  call check(all(published == 0.0_dp), 'reinitialization retained forcing values')
  call meteo%finalize()
  call restored%read('support_reordered.nc')
  call check(all(restored%ids == ids), 'rewritten restart does not follow current ID order')
  nullify(published)
  call exchange%destroy()
  call meteo%destroy()
  call meteo%destroy()
  call cache%destroy()
  call restored%destroy()
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ', message
      error stop 1
    end if
  end subroutine check
end program test_meteo_lake_restart
