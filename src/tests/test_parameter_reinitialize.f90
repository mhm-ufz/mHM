program test_parameter_reinitialize
  use mo_domain, only: domain_t
  use mo_input_container, only: input_var_dp
  use mo_grid, only: grid_t
  use mo_kind, only: i2, i4, dp
  use mo_grid_io, only: no_time
  use nml_helper, only: NML_OK
  use mo_parameter_types, only: parameter_t

  implicit none

  type(domain_t), target :: pet_domain
  type(domain_t), target :: routing_domain
  type(domain_t), target :: temperature_domain
  type(domain_t), target :: lake_domain
  type(input_var_dp) :: source
  type(grid_t), target :: source_grid
  real(dp), allocatable :: values(:)
  real(dp), allocatable :: derived_before(:)
  real(dp), pointer :: static_input(:)
  real(dp), pointer :: published_discharge(:), published_energy(:), published_temperature(:), published_lake_inflow(:)
  integer(i2), pointer :: static_fdir(:)
  character(1024) :: errmsg
  integer :: status

  call pet_domain%create(main_file="test_nml/mpr_pet_aspect_minimal.nml", cwd=".")
  status = pet_domain%exchange%config%parameters%pet_m2%set( &
    correction_factor_min=parameter_t( &
      value=0.9_dp, optimize=.true., min=0.7_dp, max=1.3_dp), &
    correction_factor_delta=parameter_t( &
      value=0.1_dp, optimize=.true., min=0.0_dp, max=0.2_dp), &
    aspect_threshold=parameter_t( &
      value=180.0_dp, optimize=.true., min=160.0_dp, max=200.0_dp), &
    errmsg=errmsg)
  call assert_true(status == NML_OK, "could not configure PET parameters through the exchange interface: " // trim(errmsg))
  call pet_domain%configure(main_file="test_nml/mpr_pet_aspect_minimal.nml")
  call pet_domain%connect()

  static_input => pet_domain%exchange%slope%data
  values = pet_domain%exchange%parameters%as_array()
  call assert_true(size(values) == 3, "PET fixture did not register the expected parameter count")
  call pet_domain%initialize()
  call assert_true(pet_domain%exchange%pet_fac_aspect%static, &
    "MPR marked an invariant parameter field as temporal")
  call assert_true(pet_domain%exchange%pet_fac_aspect%stepping == no_time, &
    "MPR did not publish an invariant parameter field with no_time stepping")
  derived_before = pet_domain%exchange%pet_fac_aspect%data

  values(1) = 1.0_dp
  call pet_domain%initialize(values)
  call assert_true( &
    maxval(abs(pet_domain%exchange%pet_fac_aspect%data - derived_before)) > 1.0e-12_dp, &
    "MPR did not rebuild its PET field after a parameter change")
  call assert_true(associated(static_input, pet_domain%exchange%slope%data), &
    "MPR repeated initialization replaced a static input")

  pet_domain%mpr%write_restart = .false.
  nullify(static_input)
  call pet_domain%finalize()
  call assert_true(associated(pet_domain%exchange%pet_fac_aspect%data), "finalize removed MPR data before destroy")
  call assert_true(allocated(pet_domain%mpr%pet%pet_fac_aspect_cache), "finalize removed MPR cache before destroy")
  call pet_domain%destroy()
  call assert_true(.not.associated(pet_domain%exchange%pet_fac_aspect%data), "destroy kept MPR publication")
  call assert_true(.not.allocated(pet_domain%mpr%pet%pet_fac_aspect_cache), "destroy kept MPR cache")
  call assert_true(.not.allocated(pet_domain%input%tgt_level0%mask), "destroy kept Input grid")
  call pet_domain%destroy()

  call routing_domain%create(main_file="test_nml/mrm_minimal.nml", cwd=".")
  call routing_domain%configure( &
    main_file="test_nml/mrm_minimal.nml", &
    para_file="test_nml/mhm_parameter_v6_routing.nml")
  call routing_domain%connect()

  published_discharge => routing_domain%exchange%discharge%data
  call assert_true(associated(published_discharge,routing_domain%mrm%flow_state%discharge), "connect did not publish flow-owned discharge")
  static_fdir => routing_domain%exchange%fdir%data
  values = routing_domain%exchange%parameters%as_array()
  call assert_true(size(values) == 1, "routing fixture did not register the expected parameter count")
  call routing_domain%initialize()
  call assert_true(associated(routing_domain%mrm%router%flow%river, routing_domain%mrm%river), "flow processor lacks river reference")
  derived_before = routing_domain%mrm%celerity
  call assert_true(allocated(routing_domain%mrm%flow_state%current_input), "mRM did not allocate flow histories")
  call assert_true(allocated(routing_domain%mrm%router%muskingum%nu1), "router did not allocate Muskingum coefficients")
  call assert_true(allocated(routing_domain%mrm%router%flow%acc_runoff), "flow processor lacks runoff scratch")

  values(1) = 2.0_dp
  call routing_domain%initialize(values)
  call assert_true(associated(routing_domain%mrm%router%flow%river, routing_domain%mrm%river), "reinitialize lost flow river reference")
  call assert_true( &
    maxval(abs(routing_domain%mrm%celerity - derived_before)) > 1.0e-12_dp, &
    "mRM did not rebuild celerity after a parameter change")
  call assert_true(associated(static_fdir, routing_domain%exchange%fdir%data), &
    "mRM repeated initialization replaced a static input")

  call assert_true(associated(published_discharge,routing_domain%mrm%flow_state%discharge), "reinitialize replaced flow discharge target")
  nullify(static_fdir,published_discharge)
  call routing_domain%finalize()
  call assert_true(associated(routing_domain%exchange%discharge%data), "finalize removed mRM discharge before destroy")
  call assert_true(allocated(routing_domain%mrm%celerity), "finalize removed mRM state before destroy")
  call routing_domain%mrm%router%destroy()
  call assert_true(.not.associated(routing_domain%mrm%router%flow%river), "processor destroy retained river reference")
  call assert_true(allocated(routing_domain%mrm%river%order%id), "processor destroy released borrowed river")
  call assert_true(allocated(routing_domain%mrm%flow_state%discharge), "router destroy released mRM-owned discharge")
  call assert_true(.not.allocated(routing_domain%mrm%router%muskingum%nu1), &
    "router destroy retained Muskingum coefficients")
  call assert_true(.not.allocated(routing_domain%mrm%router%flow%acc_runoff), &
    "router destroy retained runoff scratch")
  call routing_domain%destroy()
  call assert_true(.not.associated(routing_domain%exchange%discharge%data), "destroy kept mRM publication")
  call assert_true(.not.allocated(routing_domain%mrm%celerity), "destroy kept mRM state")
  call assert_true(.not.allocated(routing_domain%mrm%flow_state%discharge), "destroy kept router state")
  call assert_true(.not.allocated(routing_domain%mrm%level3%mask), "destroy kept mRM grid")
  call assert_true(.not.allocated(routing_domain%mrm%river%order%id), "destroy kept river order")
  call assert_true(.not.associated(routing_domain%mrm%router%river), "destroy kept router river reference")
  call assert_true(.not.associated(routing_domain%mrm%upscaler%fine_river), "destroy kept upscaler river reference")
  call assert_true(.not.allocated(routing_domain%mrm%router%scaler%id_map), "destroy kept scaler cache")
  call assert_true(.not.allocated(routing_domain%mrm%router%flow%scale_work), "destroy kept scaling work buffer")
  call assert_true(.not.associated(routing_domain%mrm%router%scaler%source_grid), "destroy kept scaler grid")
  call routing_domain%destroy()

  call temperature_domain%create(main_file="test_nml/mrm_rivertemp_meteo_minimal.nml", cwd=".")
  call temperature_domain%configure( &
    main_file="test_nml/mrm_rivertemp_meteo_minimal.nml", &
    para_file="test_nml/mhm_parameter_v6_meteo.nml")
  call temperature_domain%connect()
  published_energy => temperature_domain%exchange%energy_flux%data
  published_temperature => temperature_domain%exchange%river_temp%data
  call assert_true(associated(published_energy,temperature_domain%mrm%thermal_state%energy_flux), "connect did not publish thermal energy")
  call temperature_domain%initialize()
  call assert_true(associated(temperature_domain%mrm%router%thermal%river, temperature_domain%mrm%river), "thermal processor lacks river reference")
  call assert_true(associated(temperature_domain%exchange%energy_flux%data, &
    temperature_domain%mrm%thermal_state%energy_flux), "thermal power was not published from mRM-owned state")
  call temperature_domain%initialize()
  call assert_true(associated(temperature_domain%mrm%router%thermal%river, temperature_domain%mrm%river), "reinitialize lost thermal river reference")
  call assert_true(associated(temperature_domain%exchange%energy_flux%data, &
    temperature_domain%mrm%thermal_state%energy_flux), "reinitialize lost thermal power publication")
  call assert_true(associated(temperature_domain%exchange%river_temp%data, &
    temperature_domain%mrm%thermal_state%temperature), "reinitialize lost river temperature publication")
  call assert_true(associated(published_energy,temperature_domain%mrm%thermal_state%energy_flux), "reinitialize replaced thermal energy target")
  call assert_true(associated(published_temperature,temperature_domain%mrm%thermal_state%temperature), "reinitialize replaced thermal temperature target")
  nullify(published_energy,published_temperature)
  call temperature_domain%finalize()
  call temperature_domain%destroy()
  call assert_true(.not.associated(temperature_domain%mrm%router%thermal%river), "destroy retained thermal river reference")
  call assert_true(.not.associated(temperature_domain%exchange%energy_flux%data), &
    "destroy retained thermal power publication")

  call lake_domain%create(main_file="test_nml/mrm_lake_restart_write.nml",cwd=".")
  call lake_domain%configure(main_file="test_nml/mrm_lake_restart_write.nml", &
    para_file="test_nml/mhm_parameter_v6_routing.nml")
  lake_domain%mrm%write_restart = .false.
  lake_domain%mlm%write_restart = .false.
  call lake_domain%connect()
  published_lake_inflow => lake_domain%exchange%lake_inflow%data
  call assert_true(associated(published_lake_inflow,lake_domain%mrm%flow_state%lake_inflow), "connect did not publish flow-owned lake inflow")
  call lake_domain%initialize()
  published_lake_inflow = 7.0_dp
  call lake_domain%initialize()
  call assert_true(associated(published_lake_inflow,lake_domain%mrm%flow_state%lake_inflow), "reinitialize replaced lake inflow target")
  call assert_true(all(published_lake_inflow==0.0_dp), "fresh reinitialize retained previous lake inflow")
  call assert_true(associated(lake_domain%exchange%lake_inflow%data,lake_domain%mrm%flow_state%lake_inflow), "reinitialize lost lake exchange binding")
  nullify(published_lake_inflow)
  call lake_domain%finalize()
  call lake_domain%destroy()
  call assert_true(.not.associated(lake_domain%exchange%lake_inflow%data), "destroy retained lake inflow publication")

  allocate(source%cache(2, 2), source%mask(1, 1))
  source%grid => source_grid
  call source%destroy()
  call assert_true(.not.allocated(source%cache), "Input destroy kept dynamic cache")
  call assert_true(.not.allocated(source%mask), "Input destroy kept input mask")
  call assert_true(.not.associated(source%grid), "Input destroy kept grid reference")
  call source%destroy()

contains

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (.not.condition) then
      write(*, '(a)') trim(message)
      error stop 1
    end if
  end subroutine assert_true

end program test_parameter_reinitialize
