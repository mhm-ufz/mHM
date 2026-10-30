program test_lake_input_options
  use mo_kind, only: dp
  use mo_constants, only: nodata_dp
  use mo_domain, only: domain_t
  use mo_netcdf, only: NcDataset, NcVariable
  implicit none
  type(domain_t), target :: domain
  type(NcDataset) :: nc
  type(NcVariable) :: var
  character(:), allocatable :: name
  integer :: option
  do option = 1,2
    name = 'mlm_lake_optional'
    if (option == 2) name = 'mlm_lake_topology_input'
    nc = NcDataset(name//'.nc','a')
    if (nc%hasVariable('initial_lake_level')) then
      var = nc%getVariable('initial_lake_level')
      call var%renameVariable('unused_initial_level')
    end if
    if (option == 2) then
      if (nc%hasVariable('lake_spill_reference_surface_area')) then
        var = nc%getVariable('lake_spill_reference_surface_area')
        call var%renameVariable('unused_area')
      end if
      if (nc%hasVariable('lake_spill_reference_volume')) then
        var = nc%getVariable('lake_spill_reference_volume')
        call var%renameVariable('unused_volume')
      end if
    end if
    call nc%close()
    call domain%create(main_file='test_nml/'//name//'.nml',cwd='.')
    call domain%configure(main_file='test_nml/'//name//'.nml',para_file='test_nml/mhm_parameter_v6_lake.nml')
    domain%mlm%write_restart = .false.
    domain%meteo%lake%write_restart = .false.
    domain%mrm%write_restart = .false.
    call domain%prepare_restart()
    call domain%connect()
    call domain%initialize()
    call check(all(domain%input%lake_initial_level == nodata_dp),'optional initial level did not remain unset in Input')
    if (option == 1) then
      call check(maxval(abs(domain%exchange%lake_level%data-domain%mlm%definition%spill_reference_level-0.1_dp)) < 1.0e-12_dp, &
        'missing initial level did not use the natural-lake default')
    else
      call check(.not.domain%exchange%lake_spill_reference_surface_area%provided,'case -2 required reported bathymetry')
      call check(all(domain%mlm%definition%footprint_area > 0.0_dp),'topology-only input did not supply footprint areas')
    end if
    call domain%update()
    call domain%initialize()
    call domain%finalize()
    call domain%destroy()
  end do
contains
  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) then
      print *, 'FAIL: ',message
      error stop 1
    end if
  end subroutine check
end program test_lake_input_options
