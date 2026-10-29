program test_lake_borrowed_parameters
  use mo_kind, only: dp, i4, i8
  use mo_datetime, only: datetime
  use mo_exchange_type, only: exchange_t
  use mo_grid, only: grid_t
  use mo_points, only: points_t
  use mo_river, only: river_t
  use mo_grid_io, only: no_time
  use mo_mlm_container, only: mlm_t
  use mo_parameter_types, only: parameter_t
  implicit none
  type(exchange_t), target :: exchange
  type(grid_t), target :: grid
  type(points_t), target :: points
  type(river_t), target :: river
  type(mlm_t), target :: mlm
  type(parameter_t) :: soil(12), baseflow(10)
  character(32) :: soil_names(12), baseflow_names(10)
  integer(i8), target :: ids(1), map(2)
  integer(i4), target :: geology(2)
  real(dp), target :: level(1), area(1), volume(1), inflow(1), sand(2), clay(2)
  real(dp), allocatable :: values(:), ks(:), recession(:), borrowed(:)
  integer :: status, unit, i
  character(64) :: mode
  call get_command_argument(1,mode)
  status = exchange%config%project%set(n_geo_units=10_i4)
  call check(status == 0,'project configuration')
  status = exchange%config%processes%set(lake=1_i4,routing=2_i4,soil_moisture=3_i4,baseflow=1_i4)
  call check(status == 0,'process configuration')
  call exchange%create(cwd='.')
  call exchange%set_dims()
  exchange%cwd = exchange%root
  exchange%step_hours = 1_i4
  exchange%start_time = datetime('2001-01-01 00:00:00')
  exchange%time = exchange%start_time
  exchange%time_step_start = exchange%start_time
  call exchange%parameters%begin_configuration()
  soil = parameter_t(value=0.0_dp,optimize=.false.,min=-1000.0_dp,max=1000.0_dp)
  soil(10:12)%value = [-0.585_dp,0.0125_dp,0.0063_dp]
  baseflow = parameter_t(value=80.0_dp,optimize=.false.,min=1.0_dp,max=1000.0_dp)
  baseflow(1)%value = 10.0_dp
  do i = 1,12
    write(soil_names(i),'(a,i0)') 'coefficient',i
    if (i <= 10) write(baseflow_names(i),'(a,i0)') 'recession',i
  end do
  call exchange%parameters%add_process('soil_moisture',soil,soil_names,group='borrowed_soil')
  call exchange%parameters%add_process('baseflow',baseflow,baseflow_names,group='borrowed_baseflow')
  exchange%parameter_file = exchange%get_path('test_nml/mhm_parameter_v6_lake.nml')
  call grid%init(2_i4,1_i4,cellsize=100.0_dp)
  call points%init([50.0_dp],[50.0_dp])
  exchange%level0_lake => grid
  exchange%lake_points => points
  exchange%river_l0 => river
  ids = [101_i8]
  map = [101_i8,101_i8]
  level = 100.0_dp
  area = 1.0e4_dp
  volume = 1.0e4_dp
  inflow = 0.0_dp
  sand = [40.0_dp,70.0_dp]
  clay = [20.0_dp,10.0_dp]
  geology = [104_i4,207_i4]
  river%lake_outlet_nodes = [1_i8]
  river%upstream_area = [1.0e8_dp,1.0e8_dp]
  river%lake_map = map
  river%is_inland_sink = [.false.,.false.]
  call exchange%lake_ids%provide('Input')
  call exchange%lake_ids%prepare_data('Input',no_time)
  exchange%lake_ids%data => ids
  call exchange%lake_map%provide('Input')
  call exchange%lake_map%prepare_data('Input',no_time)
  exchange%lake_map%data => map
  call exchange%lake_spill_reference_level%provide('Input')
  call exchange%lake_spill_reference_level%prepare_data('Input',no_time)
  exchange%lake_spill_reference_level%data => level
  call exchange%lake_spill_reference_surface_area%provide('Input')
  call exchange%lake_spill_reference_surface_area%prepare_data('Input',no_time)
  exchange%lake_spill_reference_surface_area%data => area
  call exchange%lake_spill_reference_volume%provide('Input')
  call exchange%lake_spill_reference_volume%prepare_data('Input',no_time)
  exchange%lake_spill_reference_volume%data => volume
  call exchange%lake_sand%provide('Input')
  call exchange%lake_sand%prepare_data('Input',no_time)
  exchange%lake_sand%data => sand
  call exchange%lake_clay%provide('Input')
  call exchange%lake_clay%prepare_data('Input',no_time)
  exchange%lake_clay%data => clay
  call exchange%lake_geo%provide('Input')
  call exchange%lake_geo%prepare_data('Input',no_time)
  exchange%lake_geo%data => geology
  call exchange%lake_inflow%provide('mRM')
  call exchange%lake_inflow%prepare_data('mRM',1_i4)
  exchange%lake_inflow%data => inflow
  open(newunit=unit,file='lake_geology.txt',status='replace')
  write(unit,'(a)') 'nGeo_Formations 2','PARAM UNIT KARST DESCRIPTION','1 104 0 first','2 207 0 second'
  close(unit)
  status = exchange%config%lut%set(geo_lut_path=[character(256)::'lake_geology.txt'])
  call check(status == 0,'geology configuration')
  mlm%exchange => exchange
  call mlm%set_dims()
  status = mlm%config%set(n_lakes=[1_i4],bathymetry_type=reshape([1_i4],[1,1]),varying_albedo=.false.)
  call check(status == 0,'lake configuration')
  call mlm%configure()
  call exchange%parameters%seal()
  borrowed = exchange%parameters%get_process('lake')
  call check(size(borrowed) == 4,'lake duplicated borrowed PTF or recession parameters')
  call mlm%connect()
  call mlm%initialize()
  ks = mlm%definition%ks
  recession = mlm%definition%kbf
  call check(abs(recession(1)-45.0_dp) < 1.0e-12_dp,'borrowed baseflow mapping')
  mlm%write_restart = .true.
  mlm%restart_output_path = 'borrowed.nc'
  call mlm%finalize()
  if (index(mode,'frozen') == 1) then
    call exchange%destroy()
    call mlm%destroy()
    ! Recreate the same selected providers, preserving the configured registry as above.
    ! The remaining publications are supplied by the restart and the routing owner.
    call exchange%create(cwd='.')
    call exchange%set_dims()
    exchange%cwd = exchange%root
    exchange%start_time = datetime('2001-01-01 00:00:00')
    exchange%time = exchange%start_time
    exchange%step_hours = 1_i4
    call exchange%parameters%begin_configuration()
    call exchange%parameters%add_process('soil_moisture',soil,soil_names,group='borrowed_soil')
    call exchange%parameters%add_process('baseflow',baseflow,baseflow_names,group='borrowed_baseflow')
    mlm%exchange => exchange
    call mlm%set_dims()
    status = mlm%config%set(n_lakes=[1_i4],bathymetry_type=reshape([1_i4],[1,1]), &
      read_restart=[.true.],restart_input_path=[character(256)::'borrowed.nc'])
    call check(status == 0,'restart configuration')
    call mlm%configure()
    call exchange%parameters%seal()
    call mlm%prepare_restart()
    call mlm%connect()
    call exchange%lake_inflow%provide('mRM')
    call exchange%lake_inflow%prepare_data('mRM',1_i4)
    exchange%lake_inflow%data => inflow
  end if
  values = exchange%parameters%as_array()
  if (mode == 'frozen_kbf') then
    values(13) = 20.0_dp
  else
    values(10) = -0.285_dp
  end if
  call exchange%parameters%set(values)
  call mlm%initialize()
  call check(maxval(abs(mlm%definition%ks-ks)) > 1.0e-6_dp,'borrowed runtime soil override ignored')
  values(13:22) = 100.0_dp
  call exchange%parameters%set(values)
  call mlm%initialize()
  call check(abs(mlm%definition%kbf(1)-100.0_dp) < 1.0e-12_dp,'borrowed runtime baseflow override ignored')
  call exchange%destroy()
  call mlm%destroy()
contains
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    end if
  end subroutine check
end program test_lake_borrowed_parameters
