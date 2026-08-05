program test_parameter_reinitialize
  use mo_domain, only: domain_t
  use mo_kind, only: i4, dp
  use mo_grid_io, only: no_time
  use nml_helper, only: NML_OK
  use mo_parameter_types, only: parameter_t

  implicit none

  type(domain_t), target :: pet_domain
  type(domain_t), target :: routing_domain
  real(dp), allocatable :: values(:)
  real(dp), allocatable :: derived_before(:)
  real(dp), pointer :: static_input(:)
  integer(i4), pointer :: static_fdir(:)
  character(1024) :: errmsg
  integer :: status

  call pet_domain%create(main_file="test_nml/mpr_pet_aspect_minimal.nml", cwd=".")
  status = pet_domain%exchange%config%parameters%pet_m2%set( &
    correction_factor_min=parameter_t( &
      value=0.9_dp, optimize=.true., min=0.7_dp, max=1.3_dp), &
    correction_factor_max=parameter_t( &
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

  call routing_domain%create(main_file="test_nml/mrm_minimal.nml", cwd=".")
  call routing_domain%configure( &
    main_file="test_nml/mrm_minimal.nml", &
    para_file="test_nml/mhm_parameter_v6_routing.nml")
  call routing_domain%connect()

  static_fdir => routing_domain%exchange%fdir%data
  values = routing_domain%exchange%parameters%as_array()
  call assert_true(size(values) == 1, "routing fixture did not register the expected parameter count")
  call routing_domain%initialize()
  derived_before = routing_domain%mrm%river%celerity

  values(1) = 2.0_dp
  call routing_domain%initialize(values)
  call assert_true( &
    maxval(abs(routing_domain%mrm%river%celerity - derived_before)) > 1.0e-12_dp, &
    "mRM did not rebuild celerity after a parameter change")
  call assert_true(associated(static_fdir, routing_domain%exchange%fdir%data), &
    "mRM repeated initialization replaced a static input")

  nullify(static_fdir)
  call routing_domain%finalize()

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
