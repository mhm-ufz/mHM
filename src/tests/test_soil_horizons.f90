!> \file test_soil_horizons.f90
!> \brief Bounds compatibility and shared neutron-geometry regressions.
!> \date Oct 2026
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
program test_soil_horizons
  use mo_kind, only: dp, i4
  use mo_netcdf, only: NcDataset, NcDimension, NcVariable
  use mo_soil_horizons, only: read_horizon_bounds, write_horizon_metadata, compare_horizon_bounds
  use mo_mhm_constants, only: COSMIC_alpha, COSMIC_L1, COSMIC_L2, COSMIC_L4
  use mo_constants, only: PI_dp
  use mo_neutrons, only: neutron_geometry_t, DesiletsN0, COSMIC, TabularIntegralAFast
  implicit none
  type(NcDataset) :: nc
  type(NcDimension) :: dim
  type(NcVariable) :: axis, bvar
  type(neutron_geometry_t) :: work
  character(32) :: scenario
  real(dp), allocatable :: bounds(:), values(:), vertices(:), pairs(:,:)
  real(dp) :: endpoints(3), sm(3), density(3), lattice(3), l3(3), integral(10002), a, b
  real(dp) :: integral_reference, phi, attenuation, soil_mass, water_mass, expected
  integer :: n, layout, k
  logical :: variable_bottom
  real(dp), parameter :: thin_desilets(3) = [926.7530880249642_dp, 926.7530880249642_dp, 905.4058897317127_dp]
  real(dp), parameter :: thin_cosmic(3) = [804.7603257648201_dp, 1185.4589397154289_dp, 3709.1760991937017_dp]

  scenario = ''
  if (command_argument_count() > 0) call get_command_argument(1, scenario)
  do n = 1, 3
    endpoints = [100.0_dp, 300.0_dp, 1000.0_dp]
    sm = [20.0_dp, 60.0_dp, 210.0_dp]
    density = [1.2_dp, 1.4_dp, 1.6_dp]
    lattice = [1.0_dp, 2.0_dp, 3.0_dp]
    l3 = [80.0_dp, 90.0_dp, 100.0_dp]
    call work%prepare(endpoints(:n))
    call TabularIntegralAFast(integral, 20.0_dp)
    call DesiletsN0(sm(:n), endpoints(:n), density(:n), lattice(:n), 1500.0_dp, a)
    call DesiletsN0(sm(:n), endpoints(:n), density(:n), lattice(:n), 1500.0_dp, b, work)
    if (a /= b) error stop 'Desilets cached geometry changed values'
    call COSMIC(sm(:n), endpoints(:n), integral, 0.0_dp, 0.0_dp, 1500.0_dp, density(:n), lattice(:n), l3(:n), a)
    call COSMIC(sm(:n), endpoints(:n), integral, 0.0_dp, 0.0_dp, 1500.0_dp, density(:n), lattice(:n), l3(:n), b, work)
    if (a /= b) error stop 'COSMIC cached geometry changed values'
    if (n == 1) then
      ! Independent midpoint quadrature of the published COSMIC angular integral.
      soil_mass = density(1) * endpoints(1) / 20.0_dp
      water_mass = (sm(1) + lattice(1)) / 20.0_dp
      integral_reference = 0.0_dp
      do k = 1, 16384
        phi = (real(k,dp)-0.5_dp) * (PI_dp/2.0_dp) / 16384.0_dp
        attenuation = (soil_mass/l3(1) + water_mass/COSMIC_L4) / cos(phi)
        if (attenuation < 600.0_dp) integral_reference = integral_reference + exp(-attenuation)
      end do
      integral_reference = integral_reference * (PI_dp/2.0_dp) / 16384.0_dp
      expected = 1500.0_dp * exp(-soil_mass/COSMIC_L1 - water_mass/COSMIC_L2) * &
        endpoints(1)/10.0_dp * (COSMIC_alpha*density(1) + (sm(1)+lattice(1))/endpoints(1)) * &
        2.0_dp/PI_dp * integral_reference
      if (abs(a-expected) > 2.e-4_dp*expected) error stop 'COSMIC differs from independent angular quadrature'
    end if
    ! Read-only geometry is reused with changed moisture and either neutron formulation.
    sm = sm * 1.5_dp
    call COSMIC(sm(:n), endpoints(:n), integral, 2.0_dp, 3.0_dp, 1500.0_dp, density(:n), lattice(:n), l3(:n), b, work)
    call COSMIC(sm(:n), endpoints(:n), integral, 2.0_dp, 3.0_dp, 1500.0_dp, density(:n), lattice(:n), l3(:n), a)
    if (a /= b) error stop 'COSMIC retained previous moisture values'
    do layout = 0, 4
      ! New variable bounds, three old nominal formats, and current fixed bounds.
      variable_bottom = layout < 4
      vertices = [0.0_dp,endpoints(:n)]
      if (layout == 0) vertices(n+1) = vertices(n)
      nc = NcDataset('axis.nc', 'w')
      if (layout == 0 .or. layout == 4) then
        dim = nc%setCoordinate('soil', n, bounds=vertices, reference=0)
        call write_horizon_metadata(nc, 'soil', variable_bottom)
      else
        dim = nc%setCoordinate('soil', n, bounds=vertices, reference=2)
      end if
      axis = nc%getVariable('soil')
      if (layout == 2) call axis%setAttribute('bottom_boundary_type', 'domain_maximum')
      if (layout == 3) call axis%setAttribute('bottom_boundary_type', 'legacy_nominal')
      if (scenario == 'convention') call axis%setAttribute('bottom_boundary_type', 'fixed')
      if (n == 2) then
        pairs = reshape([0.0_dp,100.0_dp,100.0_dp,300.0_dp], [2,2])
        if (scenario == 'gap') pairs(1,2) = 101.0_dp
        if (scenario == 'zero-internal') pairs(:,1) = 0.0_dp
        if (scenario == 'zero-internal') pairs(1,2) = 0.0_dp
        if (scenario == 'negative-bottom') pairs(2,2) = 99.0_dp
        if (scenario == 'gap' .or. scenario == 'zero-internal' .or. scenario == 'negative-bottom') then
          bvar = nc%getVariable('soil_bnds')
          call bvar%setData(pairs)
        end if
      end if
      call nc%close()
      nc = NcDataset('axis.nc', 'r')
      call read_horizon_bounds(nc, 'soil', n, variable_bottom, bounds)
      if (variable_bottom) vertices(n+1) = vertices(n)
      if (any(bounds /= vertices)) error stop 'axis does not match normalized bounds'
      axis = nc%getVariable('soil')
      call axis%getData(values)
      if (layout == 0 .or. layout == 4) then
        if (any(values /= bounds(:n))) error stop 'coordinate does not reference layer tops'
        if (axis%hasAttribute('bottom_boundary_type')) error stop 'obsolete provenance written'
      end if
      call compare_horizon_bounds(bounds, vertices, variable_bottom)
      if (variable_bottom) call compare_horizon_bounds(bounds, [0.0_dp,endpoints(:n-1),2000.0_dp], .true.)
      if (scenario == 'internal' .and. n > 1) &
        call compare_horizon_bounds(bounds, [0.0_dp,200.0_dp,2000.0_dp], .true.)
      if (scenario == 'fixed') call compare_horizon_bounds(bounds, [0.0_dp,endpoints(:n-1),2000.0_dp], .false.)
      call nc%close()
    end do
  end do
  ! Pre-refactor numerical references: the second fractional horizon contains no 10 mm samples.
  endpoints = [11.7_dp,15.2_dp,100.5_dp]
  sm = [2.0_dp,0.5_dp,30.0_dp]
  do n = 1, 3
    call work%prepare(endpoints(:n))
    call DesiletsN0(sm(:n), endpoints(:n), density(:n), lattice(:n), 1500.0_dp, a, work)
    call COSMIC(sm(:n), endpoints(:n), integral, 2.0_dp, 3.0_dp, 1500.0_dp, density(:n), lattice(:n), l3(:n), b, work)
    if (abs(a-thin_desilets(n)) > 1.e-12_dp*thin_desilets(n)) error stop 'Desilets fractional-horizon regression'
    if (abs(b-thin_cosmic(n)) > 1.e-12_dp*thin_cosmic(n)) error stop 'COSMIC fractional-horizon regression'
  end do
  call work%destroy()
end program test_soil_horizons
