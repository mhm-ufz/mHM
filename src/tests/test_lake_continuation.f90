program test_lake_continuation
  use mo_kind, only: dp, i4
  use mo_domain, only: domain_t
  use mo_netcdf, only: NcDataset, NcVariable, NcDimension
  implicit none
  type(domain_t), target :: uninterrupted, writer, reader
  type(NcDataset) :: nc
  type(NcVariable) :: var
  type(NcDimension) :: scalar_dims(0)
  real(dp), allocatable :: values(:), bounds(:,:)
  logical, allocatable :: mask(:)
  character(256), allocatable :: names(:)
  character(64) :: mode
  character(:), allocatable :: prefix, source, topology
  integer :: lake_case, option, i, status
  call get_command_argument(1,mode)
  do lake_case = 1,2
    prefix = 'mlm_lake_case1'
    if (lake_case == 2) prefix = 'mlm_lake_case_minus2'
    if (len_trim(mode) > 0 .and. lake_case == 2) exit
    do option = 1,2
      topology = '_restart'
      if (option == 2) topology = '_topology'
      if (len_trim(mode) > 0 .and. option == 2) exit
      call writer%create(main_file='test_nml/'//prefix//topology//'_write.nml',cwd='.')
      call writer%configure(main_file='test_nml/'//prefix//topology//'_write.nml',para_file='test_nml/mhm_parameter_v6_lake.nml')
      call writer%prepare_restart()
      call writer%connect()
      call writer%initialize()
      do while(writer%exchange%time < writer%exchange%end_time)
        call writer%update()
      end do
      call writer%finalize()
      call writer%destroy()
      select case(trim(mode))
      case('unversioned')
        ! Old files are rejected before attempting to reconstruct missing definition data.
        nc = NcDataset(prefix//topology//'_mlm.nc','w')
        var = nc%setVariable('mlm_meta','i32',scalar_dims)
        call var%setData(0_i4)
        call var%setAttribute('time_stamp','1990-01-02 00:00:00')
        call nc%close()
      case('unsupported_version','case_conflict','missing_footprint','missing_ks_parameters', &
           'missing_kbf_parameters','missing_length_parameters','missing_definition')
        nc = NcDataset(prefix//topology//'_mlm.nc','a')
        select case(trim(mode))
        case('unsupported_version')
          var = nc%getVariable('mlm_meta')
          call var%setAttribute('format_version',99_i4)
        case('case_conflict')
          var = nc%getVariable('mlm_meta')
          call var%setAttribute('lake_case',-2_i4)
        case('missing_footprint')
          var = nc%getVariable('lake_footprint_area')
          call var%renameVariable('removed_footprint')
        case('missing_ks_parameters')
          var = nc%getVariable('lake_generating_ks')
          call var%renameVariable('removed_ks_parameters')
        case('missing_kbf_parameters')
          var = nc%getVariable('lake_generating_kbf')
          call var%renameVariable('removed_kbf_parameters')
        case('missing_length_parameters')
          var = nc%getVariable('lake_generating_length')
          call var%renameVariable('removed_length_parameters')
        case('missing_definition')
          var = nc%getVariable('lake_spill_reference_volume')
          call var%renameVariable('removed_definition')
        end select
        call nc%close()
      end select
      source = 'test_nml/'//prefix//topology//'_only.nml'
      if (mode == 'missing_definition') source = 'test_nml/'//prefix//'_restart_read.nml'
      call reader%create(main_file=source,cwd='.')
      if (mode == 'frozen_config_ks') then
        status = reader%exchange%config%parameters%lakes_1%from_file('test_nml/mhm_parameter_v6_lake.nml')
        call check(status == 0,'configure lake parameters')
        reader%exchange%config%parameters%lakes_1%ptf_ks_constant%value = -0.285_dp
      end if
      call reader%configure(main_file=source,para_file='test_nml/mhm_parameter_v6_lake.nml')
      if (mode == 'missing_support') reader%meteo%lake%read_restart = .false.
      if (mode == 'requested_topology') then
        reader%mlm%write_restart = .true.
        reader%mlm%restart_save_topology = .true.
      end if
      call reader%prepare_restart()
      call reader%connect()
      call reader%exchange%parameters%get_optimizer_data(values,bounds,mask,names)
      do i = 1,size(values)
        select case(trim(names(i)))
        case('lake.ptf_ks_constant')
          if (mode == 'frozen_ks') values(i) = -0.285_dp
        case('lake.spillway_length_multiplier')
          if (mode == 'frozen_length') values(i) = 2.0_dp
        end select
        if (mode == 'frozen_kbf' .and. index(names(i),'lake.baseflow_recession(') == 1) values(i)=50.0_dp
      end do
      call reader%initialize(values)
      call uninterrupted%create(main_file='test_nml/'//prefix//'_continuous.nml',cwd='.')
      call uninterrupted%configure(main_file='test_nml/'//prefix//'_continuous.nml',para_file='test_nml/mhm_parameter_v6_lake.nml')
      uninterrupted%mlm%write_restart = .false.
      uninterrupted%meteo%lake%write_restart = .false.
      uninterrupted%mrm%write_restart = .false.
      call uninterrupted%prepare_restart()
      call uninterrupted%connect()
      call uninterrupted%initialize()
      do while(uninterrupted%exchange%time < reader%exchange%start_time)
        call uninterrupted%update()
      end do
      call compare(uninterrupted,reader,lake_case)
      do while(reader%exchange%time < reader%exchange%end_time)
        call reader%update()
        call uninterrupted%update()
        call compare(uninterrupted,reader,lake_case)
      end do
      call reader%finalize()
      call reader%destroy()
      call uninterrupted%finalize()
      call uninterrupted%destroy()
    end do
  end do
contains
  subroutine compare(a,b,lake_case)
    type(domain_t), intent(in) :: a
    type(domain_t), intent(in) :: b
    integer, intent(in) :: lake_case
    ! Restored delayed coupling, lake state and routing memory must agree every hour.
    call check(maxval(abs(a%exchange%lake_outflow%data-b%exchange%lake_outflow%data)) < 1.0e-10_dp,'lake outflow continuation')
    call check(maxval(abs(a%exchange%lake_inflow%data-b%exchange%lake_inflow%data)) < 1.0e-10_dp,'delayed lake inflow continuation')
    call check(maxval(abs(a%exchange%discharge%data-b%exchange%discharge%data)) < 1.0e-10_dp,'routing continuation')
    if (lake_case == 1) then
      call check(maxval(abs(a%exchange%lake_volume%data-b%exchange%lake_volume%data)) < 1.0e-6_dp,'lake volume continuation')
      call check(maxval(abs(a%exchange%lake_gw_storage%data-b%exchange%lake_gw_storage%data)) < 1.0e-6_dp,'groundwater continuation')
      call check(maxval(abs(a%exchange%lake_level%data-b%exchange%lake_level%data)) < 1.0e-10_dp,'lake level continuation')
    end if
  end subroutine compare
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    end if
  end subroutine check
end program test_lake_continuation
