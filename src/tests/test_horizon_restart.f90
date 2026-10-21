!> \file test_horizon_restart.f90
!> \brief Restore MPR definition without soil maps/LUTs, including old v6 coordinates.
!> \date Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
program test_horizon_restart
  use mo_domain, only: domain_t
  use mo_kind, only: dp
  use mo_common_constants, only: soilHorizonsVarName
  use mo_netcdf, only: NcDataset, NcVariable
  implicit none
  type(domain_t), target :: domain
  type(NcDataset) :: nc
  type(NcVariable) :: axis, bvar, field
  real(dp), allocatable :: reference(:,:,:), bounds(:), pairs(:,:), data(:,:,:,:)
  character(32) :: scenario
  integer :: mode, legacy
  scenario = ''
  if (command_argument_count() > 0) call get_command_argument(1, scenario)
  do mode = 0, 1
    call fixture(mode, .false.)
    call start_domain()
    reference = domain%mpr%soil%sm_saturation_cache
    bounds = domain%mpr%soil%horizon_bounds
    domain%mpr%write_restart = .true.
    domain%mpr%restart_output_path = 'horizon_restart.nc'
    call domain%finalize()
    call domain%destroy()
    do legacy = 0, 3
      if (legacy > 0 .or. len_trim(scenario) > 0) then
        nc = NcDataset('horizon_restart.nc','a')
        axis = nc%getVariable(trim(soilHorizonsVarName))
        bvar = nc%getVariable(trim(soilHorizonsVarName)//'_bnds')
        call bvar%getData(pairs)
        if (legacy > 0) then
          if (mode == 0) pairs(2,2) = pairs(1,2) + 1000.0_dp
          if (legacy > 1) then
            if (mode == 1) then
              call axis%setAttribute('bottom_boundary_type','fixed')
            else if (legacy == 2) then
              call axis%setAttribute('bottom_boundary_type','domain_maximum')
            else
              call axis%setAttribute('bottom_boundary_type','legacy_nominal')
            end if
          end if
        end if
        if (scenario == 'interfaces') then
          pairs(2,1) = 150.0_dp
          pairs(1,2) = 150.0_dp
        end if
        if (scenario == 'fixed-bottom' .and. mode == 1) pairs(2,2) = pairs(2,2) + 1000.0_dp
        call bvar%setData(pairs)
        call axis%setData(pairs(2,:))
        if (scenario == 'cosmic' .and. mode == 1) then
          field = nc%getVariable('L1_COSMICL3')
          call field%getData(data)
          data = -9999.0_dp
          call field%setData(data)
        end if
        call nc%close()
      end if
      call fixture(mode, .true.)
      call start_domain()
      if (any(domain%mpr%soil%sm_saturation_cache /= reference)) error stop 'restart changed soil parameters'
      if (domain%exchange%model_horizon_bottom_variable .neqv. (mode == 0)) error stop 'wrong bottom-depth convention'
      if (any(domain%mpr%soil%horizon_bounds /= bounds)) error stop 'restart changed fixed interfaces or bottom marker'
      call domain%update()
      call domain%finalize()
      call domain%destroy()
      ! Failure scenarios for mode 1 should not edit the mode-0 file twice.
      if (len_trim(scenario) > 0) exit
    end do
  end do
contains
  !> \brief Preserve the normal lifecycle for parameter restart-only operation.
  subroutine start_domain()
    call domain%create(main_file='horizon_case.nml',cwd='.')
    call domain%configure(main_file='horizon_case.nml',para_file='test_nml/mhm_parameter_v6_hydrology.nml')
    domain%mpr%write_restart = .false.
    call domain%prepare_restart()
    call domain%connect()
    call domain%initialize()
  end subroutine start_domain

  !> \brief Build run-local cases; restart reads omit soil maps and point to absent LUTs.
  subroutine fixture(mode, restart)
    integer, intent(in) :: mode !< Soil database convention.
    logical, intent(in) :: restart !< Restore cached MPR definition.
    integer :: input, output, status
    logical :: in_input
    character(2048) :: line
    in_input = .false.
    open(newunit=input,file='test_nml/mpr_neutrons_cosmic_minimal.nml',status='old')
    open(newunit=output,file='horizon_case.nml',status='replace')
    do
      read(input,'(a)',iostat=status) line
      if (status /= 0) exit
      if (index(line,'&config_input') > 0) in_input = .true.
      if (trim(adjustl(line)) == '/') in_input = .false.
      if (restart .and. in_input) then
        if (index(line,'slope_path') > 0 .or. index(line,'aspect_path') > 0 .or. index(line,'geo_class_path') > 0 .or. &
          index(line,'soil_class_path') > 0 .or. index(line,'lai_class_path') > 0) cycle
      end if
      if (mode == 1) then
        if (index(line,'soil_db_mode(1)') > 0) line = 'soil_db_mode(1) = 1'
        if (index(line,'soil_class_path(1)') > 0) &
          line = 'soil_horizon_class_path(1) = "examples/domain_01/input/morph/soil_horizon_class.nc"'
        if (index(line,'soil_lut_path(1)') > 0) &
          line = 'soil_lut_path(1) = "examples/domain_01/input/morph/soil_classdefinition_iFlag_soilDB_1.txt"'
      end if
      if (restart) then
        if (index(line,'write_restart(1)') > 0) line = 'read_restart(1) = .true.'
        if (index(line,'restart_output_path(1)') > 0) line = 'restart_input_path(1) = "horizon_restart.nc"'
        if (index(line,'soil_lut_path(1)') > 0) line = 'soil_lut_path(1) = "absent-soil-lut.txt"'
        if (index(line,'geo_lut_path(1)') > 0) line = 'geo_lut_path(1) = "absent-geo-lut.txt"'
        if (index(line,'lai_lut_path(1)') > 0) line = 'lai_lut_path(1) = "absent-lai-lut.txt"'
        if (index(line,'land_cover_path(1)') > 0) line = 'land_cover_path(1) = "absent-land-cover.nc"'
      end if
      write(output,'(a)') trim(line)
      if (mode == 1 .and. index(line,'soil_depth(1,1)') > 0) write(output,'(a)') 'soil_depth(2,1) = 2000'
    end do
    close(input)
    close(output)
  end subroutine fixture
end program test_horizon_restart
