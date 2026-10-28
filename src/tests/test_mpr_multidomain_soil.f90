!> \file test_mpr_multidomain_soil.f90
!> \brief Regression test: every domain keeps its own MPR soil database.
!> \details The driver connects all domains before it initializes any of them. The legacy soil
!! database is module-global, so the last domain's connect used to overwrite the soil classes
!! flagged as present for the earlier domains, which then got zero soil water capacity.
!! This test builds the soil parameter fields of each domain alone and again in the driver
!! order (connect 1, connect 2, initialize 1, initialize 2) and requires identical results.
!> \authors Felix Pohl
!> \date Oct 2026
!> \changelog
!! - Felix Pohl Oct 2026: introduced the multidomain soil-database regression in e53f4fa8.
!! - Oct 2026: adapted for direct database ownership and extended lifecycle isolation coverage.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
program test_mpr_multidomain_soil
  use mo_domain, only: domain_t
  use mo_kind, only: i4, dp
  use mo_grid_scaler, only: weight_count

  implicit none

  !> \class soil_fields_t
  !> \brief Soil and neutron parameter fields of one domain.
  type :: soil_fields_t
    real(dp), allocatable :: sm_saturation(:, :, :)
    real(dp), allocatable :: sm_field_capacity(:, :, :)
    real(dp), allocatable :: wilting_point(:, :, :)
    real(dp), allocatable :: sm_exponent(:, :, :)
    real(dp), allocatable :: f_roots(:, :, :)
    real(dp), allocatable :: bulk_density(:, :, :)
    real(dp), allocatable :: lattice_water(:, :, :)
    real(dp), allocatable :: cosmic_l3(:, :, :)
  end type soil_fields_t

  character(*), parameter :: main_file = "test_nml/mpr_soil_moisture_two_domains.nml"
  character(*), parameter :: para_file = "test_nml/mhm_parameter_v6_hydrology.nml"
  type(domain_t), target :: domains(2)
  type(domain_t), target :: single
  type(soil_fields_t) :: reference(2)
  type(soil_fields_t) :: fields
  integer(i4), allocatable :: soil_ids_1(:)
  integer(i4), allocatable :: soil_ids_2(:)
  integer(i4) :: i, j, order, scenario, neutron_case
  integer(i4) :: connect_order(2)
  real(dp), allocatable :: values(:)

  do scenario = 0_i4, 2_i4
    do neutron_case = 0_i4, 2_i4
      call write_fixture(scenario, neutron_case)
      ! Standalone domains are the reference, including differing horizon counts.
      do i = 1_i4, 2_i4
        call prepare(single, i)
        call single%initialize()
        reference(i) = soil_fields(single)
        if (i == 1_i4) soil_ids_1 = present_soil_ids(single)
        if (i == 2_i4) soil_ids_2 = present_soil_ids(single)
        call single%finalize()
        call single%destroy()
      end do
      call assert_true(.not.all([(any(soil_ids_2 == soil_ids_1(i)), i = 1, size(soil_ids_1))]), &
        "fixture: domain 1 lacks an exclusive soil class")
      call assert_true(.not.all([(any(soil_ids_1 == soil_ids_2(i)), i = 1, size(soil_ids_2))]), &
        "fixture: domain 2 lacks an exclusive soil class")

      do order = 1_i4, 2_i4
        connect_order = [1_i4, 2_i4]
        if (order == 2_i4) connect_order = [2_i4, 1_i4]
        do j = 1_i4, 2_i4
          i = connect_order(j)
          call prepare(domains(i), i)
        end do
        do i = 1_i4, 2_i4
          call domains(i)%initialize()
          fields = soil_fields(domains(i))
          call assert_true(all(fields%sm_saturation > 0.0_dp), "domain has no soil water capacity")
          call compare(fields, reference(i))
        end do
        ! Parameter changes in one domain must not alter the other's input database.
        values = domains(1)%exchange%parameters%as_array()
        values(1) = values(1) * 1.1_dp
        call domains(1)%initialize(values)
        fields = soil_fields(domains(1))
        call assert_true(any(fields%sm_saturation /= reference(1)%sm_saturation), "soil parameters were not rebuilt")
        call domains(2)%initialize()
        fields = soil_fields(domains(2))
        call compare(fields, reference(2))
        call domains(1)%initialize()
        ! Explicitly restore the original vector after a changed calibration run.
        values(1) = 3.4_dp
        call domains(1)%initialize(values)
        fields = soil_fields(domains(1))
        call compare(fields, reference(1))
        call domains(1)%finalize()
        call domains(1)%destroy()
        call assert_true(.not.allocated(domains(1)%mpr%soil%db%is_present), "destroy retained soil LUT")
        call domains(1)%destroy()
        call domains(2)%initialize()
        fields = soil_fields(domains(2))
        call compare(fields, reference(2))
        call domains(2)%update()
        call domains(2)%finalize()
        call domains(2)%destroy()
      end do
    end do
  end do

  write(*, '(a)') "test_mpr_multidomain_soil: OK"

contains

  !> \brief Configure and connect a domain in the fixed lifecycle order.
  subroutine prepare(domain, id)
    type(domain_t), intent(inout), target :: domain !< Domain instance.
    integer(i4), intent(in) :: id !< Configured domain index.
    call domain%create(main_file="soil_case.nml", domain=id, cwd=".")
    call domain%configure(main_file="soil_case.nml", para_file=para_file)
    domain%mpr%write_restart = .false.
    call domain%prepare_restart()
    call domain%connect()
    call assert_true(domain%mpr%upscaler%weight_mode == weight_count, "MPR scaler must use count weights")
    if (domain%mpr%soil%db%mode == 0) then
      call assert_true(domain%exchange%model_horizon_bottom_variable, 'missing variable-bottom flag')
      call assert_true(domain%mpr%soil%horizon_bounds(domain%exchange%n_layers+1) == &
        domain%mpr%soil%horizon_bounds(domain%exchange%n_layers), 'variable bottom does not repeat its top')
    else
      call assert_true(.not.domain%exchange%model_horizon_bottom_variable, 'missing fixed-bottom flag')
    end if
  end subroutine prepare

  !> \brief Derive small scenario namelists from the two-domain regression fixture.
  subroutine write_fixture(scenario, neutron_case)
    integer(i4), intent(in) :: scenario !< Conventional, mixed, or single-layer soil layouts.
    integer(i4), intent(in) :: neutron_case !< Disabled, Desilets, or COSMIC neutrons.
    integer :: input_unit, output_unit, status
    character(2048) :: line
    open(newunit=input_unit, file=main_file, status="old", action="read")
    open(newunit=output_unit, file="soil_case.nml", status="replace", action="write")
    do
      read(input_unit, '(a)', iostat=status) line
      if (status /= 0) exit
      if (index(line, "n_layers(2)") > 0) line = "  n_layers(2) = 3"
      if (scenario == 2_i4 .and. index(line, "n_layers(1)") > 0) line = "  n_layers(1) = 1"
      if (scenario == 1_i4) then
        if (index(line, "soil_db_mode(1)") > 0) line = "  soil_db_mode(1) = 1"
        if (index(line, "soil_class_path(1)") > 0) &
          line = '  soil_horizon_class_path(1) = "examples/domain_01/input/morph/soil_horizon_class.nc"'
        if (index(line, "soil_lut_path(1)") > 0) &
          line = '  soil_lut_path(1) = "examples/domain_01/input/morph/soil_classdefinition_iFlag_soilDB_1.txt"'
      end if
      write(output_unit, '(a)') trim(line)
      if (index(line, "soil_moisture =") > 0) write(output_unit, '(a,i0)') "  neutrons = ", neutron_case
      if (index(line, "soil_depth(1,2)") > 0) write(output_unit, '(a)') "  soil_depth(2,2) = 500"
      if (scenario == 1_i4 .and. index(line, "soil_depth(1,1)") > 0) &
        write(output_unit, '(a)') "  soil_depth(2,1) = 2000"
    end do
    close(input_unit)
    close(output_unit)
  end subroutine write_fixture

  !> \brief Compare domain-owned fields against standalone references.
  subroutine compare(actual, expected)
    type(soil_fields_t), intent(in) :: actual !< Current fields.
    type(soil_fields_t), intent(in) :: expected !< Standalone reference fields.
    call assert_same(actual%sm_saturation, expected%sm_saturation, "sm_saturation")
    call assert_same(actual%sm_field_capacity, expected%sm_field_capacity, "sm_field_capacity")
    call assert_same(actual%wilting_point, expected%wilting_point, "wilting_point")
    call assert_same(actual%sm_exponent, expected%sm_exponent, "sm_exponent")
    call assert_same(actual%f_roots, expected%f_roots, "f_roots")
    call assert_same(actual%bulk_density, expected%bulk_density, "bulk_density")
    call assert_same(actual%lattice_water, expected%lattice_water, "lattice_water")
    call assert_same(actual%cosmic_l3, expected%cosmic_l3, "cosmic_l3")
  end subroutine compare

  !> \brief Copy the cached soil parameter fields of a domain.
  function soil_fields(domain) result(fields)
    type(domain_t), intent(in) :: domain !< Domain supplying the fields.
    type(soil_fields_t) :: fields

    fields%sm_saturation = domain%mpr%soil%sm_saturation_cache
    fields%sm_field_capacity = domain%mpr%soil%sm_field_capacity_cache
    fields%wilting_point = domain%mpr%soil%wilting_point_cache
    fields%sm_exponent = domain%mpr%soil%sm_exponent_cache
    fields%f_roots = domain%mpr%soil%f_roots_cache
    if (allocated(domain%mpr%neutron%bulk_density_cache)) then
      fields%bulk_density = domain%mpr%neutron%bulk_density_cache
    else
      allocate(fields%bulk_density(0, 0, 0))
    end if
    if (allocated(domain%mpr%neutron%lattice_water_cache)) then
      fields%lattice_water = domain%mpr%neutron%lattice_water_cache
    else
      allocate(fields%lattice_water(0, 0, 0))
    end if
    if (allocated(domain%mpr%neutron%cosmic_l3_cache)) then
      fields%cosmic_l3 = domain%mpr%neutron%cosmic_l3_cache
    else
      allocate(fields%cosmic_l3(0, 0, 0))
    end if
  end function soil_fields

  !> \brief Distinct soil class ids of the first soil layer of a domain.
  function present_soil_ids(domain) result(ids)
    type(domain_t), intent(in) :: domain !< Domain supplying the fields.
    integer(i4), allocatable :: ids(:)
    integer(i4) :: k
    integer(i4) :: id

    allocate(ids(0))
    do k = 1_i4, size(domain%exchange%soil_id%data, 1)
      id = domain%exchange%soil_id%data(k, 1)
      if (id < 1_i4) cycle
      if (.not.any(ids == id)) ids = [ids, id]
    end do
  end function present_soil_ids

  !> \brief Require exactly reproducible parameter fields.
  subroutine assert_same(actual, expected, name)
    real(dp), intent(in) :: actual(:, :, :) !< Computed fields.
    real(dp), intent(in) :: expected(:, :, :) !< Reference fields.
    character(*), intent(in) :: name !< Field name.

    call assert_true(all(shape(actual) == shape(expected)), name // " shape differs from the single-domain run")
    call assert_true(all(actual == expected), name // " differs from the single-domain run")
  end subroutine assert_same

  !> \brief Terminate on a failed regression invariant.
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition !< Required invariant.
    character(*), intent(in) :: message !< Failure diagnostic.

    if (.not.condition) then
      write(*, '(a)') trim(message)
      error stop 1
    end if
  end subroutine assert_true

end program test_mpr_multidomain_soil
