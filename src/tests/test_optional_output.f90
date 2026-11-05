!> \file test_optional_output.f90
!> \brief Compare omitted and empty output groups, including repeated reads and failures.
!> \version 0.1
!> \date November 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
program test_optional_output
  use mo_kind, only: i4
  use nml_helper, only: NML_OK, NML_ERR_FILE_NOT_FOUND, NML_ERR_READ, NML_ERR_NML_NOT_FOUND
  use nml_baseflow_1, only: nml_baseflow_1_t
  use nml_output_mhm, only: nml_output_mhm_t
  use nml_output_mlm, only: nml_output_mlm_t
  use nml_output_mrm, only: nml_output_mrm_t
  implicit none
  integer :: unit

  open(newunit=unit, file='omitted.nml', status='replace')
  write(unit, '(a)') '! No output groups.', '&unrelated', '/'
  close(unit)
  open(newunit=unit, file='empty.nml', status='replace')
  write(unit, '(a)') '&output_mhm', '/', '&output_mlm', '/', '&output_mrm', '/'
  close(unit)
  open(newunit=unit, file='malformed.nml', status='replace')
  write(unit, '(a)') '&output_mhm', 'unknown_option=1', '/', &
    '&output_mlm', 'unknown_option=1', '/', '&output_mrm', 'unknown_option=1', '/'
  close(unit)
  call check_mhm()
  call check_mlm()
  call check_mrm()
  call check_required()

contains

  !> \brief Stop with a useful diagnostic when a reader contract is violated.
  subroutine check(condition, message)
    logical, intent(in) :: condition !< Expected condition.
    character(*), intent(in) :: message !< Failure diagnostic.
    if (.not.condition) then
      print *, message
      error stop 1
    end if
  end subroutine check

  !> \brief A group with required parameter values must still be present.
  subroutine check_required()
    type(nml_baseflow_1_t) :: config
    integer :: status
    character(1024) :: errmsg
    status = config%set_dims(n_geo_units=2, errmsg=errmsg)
    call check(status == NML_OK, 'could not size required parameter group: '//trim(errmsg))
    status = config%from_file('omitted.nml', errmsg=errmsg)
    call check(status == NML_ERR_NML_NOT_FOUND, 'missing required parameter group did not fail')
    call check(.not.config%is_configured, 'missing required parameter group was configured')
    call check(len_trim(errmsg) > 0, 'missing required parameter group lacked a diagnostic')
  end subroutine check_required

  !> \brief Check mhm output defaults and reset/error behavior through the public reader.
  subroutine check_mhm()
    type(nml_output_mhm_t) :: config
    integer :: status, mode
    character(1024) :: errmsg
    character(11) :: file

    do mode = 1, 2
      file = 'omitted.nml'
      if (mode == 2) file = 'empty.nml'
      ! Prime every switch so a later read must reset the entire configuration.
      status = config%set(output_frequency=72_i4, output_deflate_level=1_i4, &
        output_time_reference=0_i4, output_double_precision=.true., &
        out_interception=.true., &
        out_snowpack=.true., &
        out_swc=.true., &
        out_sm=.true., &
        out_sm_all=.true., &
        out_sealedstw=.true., &
        out_unsatstw=.true., &
        out_satstw=.true., &
        out_pet=.true., &
        out_aet_all=.true., &
        out_q=.true., &
        out_qd=.true., &
        out_qif=.true., &
        out_qis=.true., &
        out_qb=.true., &
        out_recharge=.true., &
        out_soil_infil=.true., &
        out_neutrons=.true., &
        out_aet_layer=.true., &
        out_preeffect=.true., &
        out_qsm=.true., errmsg=errmsg)
      call check(status == NML_OK, 'could not prime mhm output: '//trim(errmsg))
      status = config%from_file(trim(file), errmsg=errmsg)
      call check(status == NML_OK, 'mhm optional read failed: '//trim(errmsg))
      call check(config%is_configured, 'mhm optional read did not configure defaults')
      call check(len_trim(errmsg) == 0, 'mhm successful read retained an error message')
      call check(config%output_frequency == -2_i4, 'mhm frequency default changed')
      call check(config%output_deflate_level == 6_i4, 'mhm compression default changed')
      call check(config%output_time_reference == 2_i4, 'mhm timestamp default changed')
      call check(.not.config%output_double_precision, 'mhm precision default changed')
      call check(.not.any([ &
        config%out_interception, &
        config%out_snowpack, &
        config%out_swc, &
        config%out_sm, &
        config%out_sm_all, &
        config%out_sealedstw, &
        config%out_unsatstw, &
        config%out_satstw, &
        config%out_pet, &
        config%out_aet_all, &
        config%out_q, &
        config%out_qd, &
        config%out_qif, &
        config%out_qis, &
        config%out_qb, &
        config%out_recharge, &
        config%out_soil_infil, &
        config%out_neutrons, &
        config%out_aet_layer, &
        config%out_preeffect, &
        config%out_qsm]), 'mhm switches did not reset to defaults')
      status = config%is_valid(errmsg=errmsg)
      call check(status == NML_OK, 'mhm optional defaults invalid: '//trim(errmsg))
    end do
    status = config%from_file('missing-output-file.nml', errmsg=errmsg)
    call check(status == NML_ERR_FILE_NOT_FOUND, 'mhm missing file did not fail')
    call check(len_trim(errmsg) > 0, 'mhm missing file lacked a diagnostic')
    status = config%from_file('malformed.nml', errmsg=errmsg)
    call check(status == NML_ERR_READ, 'mhm malformed group did not fail')
    call check(len_trim(errmsg) > 0, 'mhm malformed group lacked a diagnostic')
  end subroutine check_mhm

  !> \brief Check mlm output defaults and reset/error behavior through the public reader.
  subroutine check_mlm()
    type(nml_output_mlm_t) :: config
    integer :: status, mode
    character(1024) :: errmsg
    character(11) :: file

    do mode = 1, 2
      file = 'omitted.nml'
      if (mode == 2) file = 'empty.nml'
      ! Prime every switch so a later read must reset the entire configuration.
      status = config%set(output_frequency=72_i4, output_deflate_level=1_i4, &
        output_time_reference=0_i4, output_double_precision=.true., &
        output_timeseries=.true., &
        out_lake_level=.true., &
        out_lake_surface_area=.true., &
        out_lake_volume=.true., &
        out_lake_gw_storage=.true., &
        out_lake_inflow=.true., &
        out_lake_precipitation=.true., &
        out_lake_evaporation=.true., &
        out_lake_seepage=.true., &
        out_lake_baseflow=.true., &
        out_lake_spill=.true., &
        out_lake_outflow=.true., errmsg=errmsg)
      call check(status == NML_OK, 'could not prime mlm output: '//trim(errmsg))
      status = config%from_file(trim(file), errmsg=errmsg)
      call check(status == NML_OK, 'mlm optional read failed: '//trim(errmsg))
      call check(config%is_configured, 'mlm optional read did not configure defaults')
      call check(len_trim(errmsg) == 0, 'mlm successful read retained an error message')
      call check(config%output_frequency == -1_i4, 'mlm frequency default changed')
      call check(config%output_deflate_level == 6_i4, 'mlm compression default changed')
      call check(config%output_time_reference == 2_i4, 'mlm timestamp default changed')
      call check(.not.config%output_double_precision, 'mlm precision default changed')
      call check(.not.any([ &
        config%output_timeseries, &
        config%out_lake_level, &
        config%out_lake_surface_area, &
        config%out_lake_volume, &
        config%out_lake_gw_storage, &
        config%out_lake_inflow, &
        config%out_lake_precipitation, &
        config%out_lake_evaporation, &
        config%out_lake_seepage, &
        config%out_lake_baseflow, &
        config%out_lake_spill, &
        config%out_lake_outflow]), 'mlm switches did not reset to defaults')
      status = config%is_valid(errmsg=errmsg)
      call check(status == NML_OK, 'mlm optional defaults invalid: '//trim(errmsg))
    end do
    status = config%from_file('missing-output-file.nml', errmsg=errmsg)
    call check(status == NML_ERR_FILE_NOT_FOUND, 'mlm missing file did not fail')
    call check(len_trim(errmsg) > 0, 'mlm missing file lacked a diagnostic')
    status = config%from_file('malformed.nml', errmsg=errmsg)
    call check(status == NML_ERR_READ, 'mlm malformed group did not fail')
    call check(len_trim(errmsg) > 0, 'mlm malformed group lacked a diagnostic')
  end subroutine check_mlm

  !> \brief Check mrm output defaults and reset/error behavior through the public reader.
  subroutine check_mrm()
    type(nml_output_mrm_t) :: config
    integer :: status, mode
    character(1024) :: errmsg
    character(11) :: file

    do mode = 1, 2
      file = 'omitted.nml'
      if (mode == 2) file = 'empty.nml'
      ! Prime every switch so a later read must reset the entire configuration.
      status = config%set(output_frequency=72_i4, output_deflate_level=1_i4, &
        output_time_reference=0_i4, output_double_precision=.true., &
        out_facc=.true., &
        out_upstream_area=.true., &
        out_reach_surface_area=.true., &
        out_qrouted=.true., &
        out_rivtemp=.true., &
        out_erouted=.true., errmsg=errmsg)
      call check(status == NML_OK, 'could not prime mrm output: '//trim(errmsg))
      status = config%from_file(trim(file), errmsg=errmsg)
      call check(status == NML_OK, 'mrm optional read failed: '//trim(errmsg))
      call check(config%is_configured, 'mrm optional read did not configure defaults')
      call check(len_trim(errmsg) == 0, 'mrm successful read retained an error message')
      call check(config%output_frequency == -1_i4, 'mrm frequency default changed')
      call check(config%output_deflate_level == 6_i4, 'mrm compression default changed')
      call check(config%output_time_reference == 2_i4, 'mrm timestamp default changed')
      call check(.not.config%output_double_precision, 'mrm precision default changed')
      call check(.not.any([ &
        config%out_facc, &
        config%out_upstream_area, &
        config%out_reach_surface_area, &
        config%out_qrouted, &
        config%out_rivtemp, &
        config%out_erouted]), 'mrm switches did not reset to defaults')
      status = config%is_valid(errmsg=errmsg)
      call check(status == NML_OK, 'mrm optional defaults invalid: '//trim(errmsg))
    end do
    status = config%from_file('missing-output-file.nml', errmsg=errmsg)
    call check(status == NML_ERR_FILE_NOT_FOUND, 'mrm missing file did not fail')
    call check(len_trim(errmsg) > 0, 'mrm missing file lacked a diagnostic')
    status = config%from_file('malformed.nml', errmsg=errmsg)
    call check(status == NML_ERR_READ, 'mrm malformed group did not fail')
    call check(len_trim(errmsg) > 0, 'mrm malformed group lacked a diagnostic')
  end subroutine check_mrm

end program test_optional_output
