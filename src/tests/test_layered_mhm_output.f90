!> \file test_layered_mhm_output.f90
!> \brief Layered/ordinary output values, precision, aggregation, and singleton axes.
!> \date Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
program test_layered_mhm_output
  use mo_domain, only: domain_t
  use mo_kind, only: dp, i4
  use mo_netcdf, only: NcDataset, NcVariable
  implicit none
  type(domain_t), target :: domain
  type(NcDataset) :: nc
  type(NcVariable) :: field
  real(dp), allocatable :: expected(:,:,:,:), ordinary(:,:), packed(:), actual(:,:,:,:), flat(:,:,:), tops(:), pairs(:,:)
  character(16), parameter :: names(4) = [character(16) :: 'SWC','SM','soil_infil','aET_soil']
  character(256) :: comment
  integer :: nh, mode, layout, precision, step, i, j, day, h, k, idx, n, unit, input, status
  character(2048) :: line

  do layout = 1, 3
    nh = min(layout, 2)
    mode = 0
    if (layout == 3) mode = 1
    do precision = 1, 2
      do idx = 1, 2
        step = 1
        if (idx == 2) step = 24
        open(newunit=input, file='test_nml/mhm_output_minimal.nml', status='old')
        open(newunit=unit, file='layered_case.nml', status='replace')
        do
          read(input,'(a)',iostat=status) line
          if (status /= 0) exit
          if (index(line,'n_layers(1)') > 0) write(line,'(a,i0)') 'n_layers(1) = ',nh
          if (index(line,'neutrons =') > 0) line = 'neutrons = 0'
          if (mode == 1) then
            if (index(line,'soil_db_mode(1)') > 0) line = 'soil_db_mode(1) = 1'
            if (index(line,'soil_class_path(1)') > 0) &
              line = 'soil_horizon_class_path(1) = "examples/domain_01/input/morph/soil_horizon_class.nc"'
            if (index(line,'soil_lut_path(1)') > 0) &
              line = 'soil_lut_path(1) = "examples/domain_01/input/morph/soil_classdefinition_iFlag_soilDB_1.txt"'
          end if
          write(unit,'(a)') trim(line)
          if (index(line,'&config_time') > 0) write(unit,'(a,i0)') 'time_step(1) = ',step
          if (mode == 1 .and. index(line,'soil_depth(1,1)') > 0) write(unit,'(a)') 'soil_depth(2,1) = 2000'
        end do
        close(input)
        close(unit)
        open(newunit=unit, file='layered_out.nml', status='replace')
        write(unit,'(a)') '&output_mhm', 'output_frequency = -1', 'out_swc = .true.', &
          'out_sm = .true.', 'out_soil_infil = .true.', 'out_aet_layer = .true.', &
          'out_sm_all = .true.', 'out_aet_all = .true.'
        write(unit,'(a,l1)') 'output_double_precision = ',precision == 2
        write(unit,'(a)') '/'
        close(unit)
        call domain%create(main_file='layered_case.nml', cwd='.')
        call domain%configure(main_file='layered_case.nml', para_file='test_nml/mhm_parameter_v6_hydrology.nml', &
          out_file='layered_out.nml')
        domain%mpr%write_restart = .false.
        call domain%prepare_restart()
        call domain%connect()
        call domain%initialize()
        n = int(domain%exchange%level1_land%ncells)
        allocate(expected(n,nh,2,4),ordinary(n,2),packed(n))
        expected = 0.0_dp
        ordinary = 0.0_dp
        do i = 1, 48/step
          call domain%update()
          day = (i-1)/(24/step)+1
          expected(:,:,day,1) = expected(:,:,day,1) + domain%exchange%soil_moisture%data / real(24/step,dp)
          expected(:,:,day,2) = expected(:,:,day,2) + &
            domain%exchange%soil_moisture%data / domain%exchange%sm_saturation%data / real(24/step,dp)
          do h = 1, nh
            expected(:,h,day,3) = expected(:,h,day,3) + &
              domain%exchange%infiltration%data(:,h) * (1.0_dp-domain%exchange%f_sealed%data)
            expected(:,h,day,4) = expected(:,h,day,4) + &
              domain%exchange%aet_soil%data(:,h) * (1.0_dp-domain%exchange%f_sealed%data)
          end do
          ordinary(:,day) = ordinary(:,day) + sum(domain%exchange%soil_moisture%data,dim=2) / &
            sum(domain%exchange%sm_saturation%data,dim=2) / real(24/step,dp)
        end do
        call domain%finalize()
        nc = NcDataset('mhm_output_minimal.nc','r')
        if (nc%hasVariable('SWC_L01')) error stop 'legacy numbered output retained'
        do j = 1, 4
          field = nc%getVariable(trim(names(j)))
          call field%getData(actual)
          if (size(actual,3) /= nh .or. size(actual,4) /= 2) error stop 'wrong layered output shape'
          do day = 1, 2
            do h = 1, nh
              call domain%exchange%level1_land%pack_into(actual(:,:,h,day), packed)
              do k = 1, n
                if (abs(packed(k)-expected(k,h,day,j)) > 1.e-6_dp * max(1.0_dp,abs(expected(k,h,day,j)))) &
                  error stop 'layered output value or aggregation changed'
              end do
            end do
          end do
        end do
        field = nc%getVariable('SM_Lall')
        call field%getData(flat)
        do day = 1, 2
          call domain%exchange%level1_land%pack_into(flat(:,:,day),packed)
          if (any(abs(packed-ordinary(:,day)) > 1.e-6_dp)) error stop 'ordinary output changed'
        end do
        field = nc%getVariable('z')
        call field%getData(tops)
        if (any(tops /= domain%mpr%soil%horizon_bounds(:nh))) error stop 'layer tops not written'
        if (field%hasAttribute('bottom_boundary_type')) error stop 'obsolete provenance attribute retained'
        if (mode == 0) then
          call field%getAttribute('comment',comment)
          if (index(comment,'mHM convention: coincident bounds') == 0) error stop 'missing variable-bottom explanation'
        else
          if (field%hasAttribute('comment')) error stop 'fixed bounds described as variable'
        end if
        field = nc%getVariable('z_bnds')
        call field%getData(pairs)
        if (any(pairs(1,:) /= tops)) error stop 'output bounds do not match layer tops'
        if (any(pairs(2,:) /= domain%mpr%soil%horizon_bounds(2:))) error stop 'output bounds changed'
        if (mode == 0) then
          if (pairs(1,nh) /= pairs(2,nh)) error stop 'variable bottom has an invented thickness'
        else
          if (any(pairs(2,:) <= pairs(1,:))) error stop 'fixed horizon has nonpositive thickness'
        end if
        call nc%close()
        call domain%destroy()
        deallocate(expected,ordinary,packed)
      end do
    end do
  end do
end program test_layered_mhm_output
