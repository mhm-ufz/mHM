program test_mlm_parameters
  use mo_kind, only: dp
  use mo_domain, only: domain_t
  implicit none
  type(domain_t), target :: domain, independent
  real(dp), allocatable :: values(:), bounds(:,:), old_ks(:), old_kbf(:), old_length(:)
  logical, allocatable :: mask(:)
  character(256), allocatable :: names(:)
  real(dp), pointer :: level(:), forcing(:)
  integer :: i, status
  call domain%create(main_file='test_nml/mlm_lake_case1_pt.nml',cwd='.')
  status = domain%exchange%config%parameters%lakes_1%from_file('test_nml/mhm_parameter_v6_lake.nml')
  call check(status == 0, 'interface parameter setup')
  domain%exchange%config%parameters%lakes_1%spill_coefficient_of_discharge%value = 2.2_dp
  call domain%configure(main_file='test_nml/mlm_lake_case1_pt.nml',para_file='test_nml/mhm_parameter_v6_lake.nml')
  call domain%connect()
  call domain%exchange%parameters%get_optimizer_data(values,bounds,mask,names)
  call domain%initialize(values)
  call check(abs(domain%mlm%runtime%cd-2.2_dp) < 1.0e-12_dp,'parameter file overwrote interface configuration')
  old_ks = domain%mlm%definition%ks
  old_kbf = domain%mlm%definition%kbf
  old_length = domain%mlm%definition%spillway_length
  level => domain%exchange%lake_level%data
  forcing => domain%exchange%lake_temp%data
  ! Change the paths after connect: initialization must consume cached static inputs.
  domain%exchange%config%lut%geo_lut_path(1) = 'must_not_be_read.txt'
  domain%exchange%config%lut%soil_lut_path(1) = 'must_not_be_read.txt'
  do i = 1,size(values)
    if (index(names(i), "lake.baseflow_recession(") == 1) values(i)=50.0_dp
    select case(trim(names(i)))
    case('lake.spill_coefficient_of_discharge'); values(i)=1.8_dp
    case('lake.spillway_length_multiplier'); values(i)=2.0_dp
    case('lake.b'); values(i)=250.0_dp
    case('lake.ptf_ks_constant'); values(i)=-0.285_dp
    case('lake.baseflow_recession(1)'); values(i)=50.0_dp
    end select
  end do
  call domain%initialize(values)
  call check(associated(level,domain%exchange%lake_level%data),'reinitialization replaced lake state publication')
  call check(associated(forcing,domain%exchange%lake_temp%data),'reinitialization replaced forcing publication')
  call check(maxval(abs(domain%mlm%definition%ks-old_ks)) > 1.0e-6_dp,'runtime Cosby override ignored')
  call check(maxval(abs(domain%mlm%definition%kbf-old_kbf)) > 1.0e-6_dp,'runtime recession override ignored')
  call check(maxval(abs(domain%mlm%definition%spillway_length-2.0_dp*old_length)) < 1.0e-10_dp,'runtime length override ignored')
  call check(abs(domain%mlm%runtime%b-250.0_dp) < 1.0e-12_dp,'runtime thickness override ignored')
  call domain%update()
  call domain%initialize(values)
  call check(all(domain%mlm%state%spill == 0.0_dp),'reinitialization retained spill diagnostics')
  call check(all(domain%mlm%state%evaporation == 0.0_dp),'reinitialization retained evaporation diagnostics')
  call check(all(forcing == 0.0_dp),'reinitialization retained forcing')
  call independent%create(main_file='test_nml/mlm_lake_case1_pt.nml',cwd='.')
  call independent%configure(main_file='test_nml/mlm_lake_case1_pt.nml',para_file='test_nml/mhm_parameter_v6_lake.nml')
  call independent%connect()
  call independent%initialize()
  call check(maxval(abs(independent%mlm%definition%ks-old_ks)) < 1.0e-12_dp,'lake parameters leaked between domains')
  call check(.not.associated(independent%exchange%lake_level%data,level),'domains share mutable state')
  call independent%finalize()
  call independent%destroy()
  nullify(level,forcing)
  call domain%finalize()
  call domain%destroy()
  call check(.not.allocated(domain%meteo%lake%temp),'destroy retained lake temperature')
  call domain%create(main_file='test_nml/mlm_lake_case1_pt.nml',cwd='.')
  call domain%configure(main_file='test_nml/mlm_lake_case1_pt.nml',para_file='test_nml/mhm_parameter_v6_lake.nml')
  call domain%connect()
  call domain%initialize()
  call domain%finalize()
  call domain%destroy()
contains
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    end if
  end subroutine check
end program test_mlm_parameters
