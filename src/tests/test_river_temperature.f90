!> Independent heat-balance checks for the open-water routing kernel.
program test_river_temperature
  use mo_kind, only: dp, i4
  use mo_constants, only: T0_dp, sigma_dp, Psychro_dp
  use mo_pet, only: slope_satpressure
  use mo_river_temperature, only: water_heat_capacity_volume, river_heat_flux, solve_river_temperature
  implicit none
  real(dp) :: t, source, q, area, c2, d, ta, sw, lw, albedo, pt, emissivity, kh, f, h0
  real(dp) :: seed, target, eta, a4
  integer(i4) :: status, case_id
  logical :: fallback

  if (abs(water_heat_capacity_volume - 4.19e6_dp) > 1.0e-6_dp) error stop 'water volumetric heat capacity'
  q = 2.0_dp
  call solve_river_temperature(q, water_heat_capacity_volume*q*280.0_dp, 0.0_dp, 0.0_dp, &
    280.0_dp, 285.0_dp, 0.0_dp, 0.0_dp, 0.08_dp, 1.26_dp, 0.98_dp, 20.0_dp, &
    t, source, status, fallback)
  if (status /= 0_i4 .or. abs(t - 280.0_dp) > 1.0e-6_dp .or. abs(source) > 1.0e-8_dp) &
    error stop 'linear heat balance'

  call solve_river_temperature(0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 280.0_dp, T0_dp, &
    0.0_dp, 0.0_dp, 0.08_dp, 1.26_dp, 0.98_dp, 20.0_dp, t, source, status, fallback)
  if (status /= 1_i4 .or. fallback) error stop 'reject nonpositive solver discharge'
  call solve_river_temperature(q, 0.0_dp, 0.0_dp, -1.0_dp, 280.0_dp, T0_dp, &
    0.0_dp, 0.0_dp, 0.08_dp, 1.26_dp, 0.98_dp, 20.0_dp, t, source, status, fallback)
  if (status /= 1_i4 .or. fallback) error stop 'reject negative reach area'

  area = 1000.0_dp
  c2 = 0.5_dp
  d = c2 * area * 0.5_dp
  ta = 290.0_dp
  albedo = 0.08_dp
  emissivity = 0.98_dp
  kh = 20.0_dp
  do case_id = 1_i4, 2_i4
    if (case_id == 1_i4) then
      sw = 500.0_dp
      lw = 300.0_dp
      pt = 2.0_dp ! permits a negative quartic coefficient on the sunlit branch
    else
      sw = 0.0_dp
      lw = 100.0_dp
      pt = 1.26_dp
    end if
    seed = 260.0_dp
    target = 280.0_dp
    f = 0.5_dp * (river_heat_flux(seed, ta, sw, lw, albedo, pt, emissivity, kh) + &
      river_heat_flux(target, ta, sw, lw, albedo, pt, emissivity, kh))
    h0 = water_heat_capacity_volume*q*target - c2*area*f
    call solve_river_temperature(q, h0, d, area, seed, ta, sw, lw, &
      albedo, pt, emissivity, kh, t, source, status, fallback)
    if (status /= 0_i4 .or. abs(t - 280.0_dp) > 1.0e-6_dp) error stop 'quartic heat balance'
    if (abs(water_heat_capacity_volume*q*t - h0 - c2*source) > 1.0e-4_dp) &
      error stop 'physical energy balance in watts'
    if (abs(source - area*f) > 1.0e-6_dp) error stop 'surface source in watts'
  end do

  ! Exact branch boundary and subfreezing results exercise the original clamped heat balance.
  pt = 1.26_dp
  sw = 500.0_dp
  lw = 300.0_dp
  target = (((1.0_dp-albedo)*sw+lw)/(emissivity*sigma_dp))**0.25_dp
  call check_known_root(target, 250.0_dp)
  sw = 0.0_dp
  lw = 100.0_dp
  call check_known_root(260.0_dp, 290.0_dp)

  ! A sunlit negative quartic has a maximum at 320 K; splitting there must retain both roots.
  ta = 310.0_dp
  sw = 2000.0_dp
  lw = 300.0_dp
  pt = 2.0_dp
  kh = 0.0_dp
  eta = pt*slope_satpressure(ta-T0_dp)/(Psychro_dp+slope_satpressure(ta-T0_dp))
  a4 = d*(1.0_dp-eta)*emissivity*sigma_dp
  q = -4.0_dp*a4*320.0_dp**3/water_heat_capacity_volume
  seed = 320.0_dp
  target = 280.0_dp
  h0 = water_heat_capacity_volume*q*target - d * &
    (river_heat_flux(seed,ta,sw,lw,albedo,pt,emissivity,kh) + &
     river_heat_flux(target,ta,sw,lw,albedo,pt,emissivity,kh))
  call solve_river_temperature(q,h0,d,area,seed,ta,sw,lw,albedo,pt,emissivity,kh,t,source,status,fallback)
  if (status /= 2_i4) error stop 'multiple admissible roots'
  if (.not.fallback) error stop 'stationary seed must trigger bracket safeguard'
  call check_known_root(320.0_dp, 320.0_dp)

  call solve_river_temperature(q, -1.0e12_dp, 0.0_dp, 0.0_dp, 280.0_dp, T0_dp, &
    0.0_dp, 0.0_dp, albedo, pt, emissivity, kh, t, source, status, fallback)
  if (status /= 1_i4) error stop 'reject root outside numerical interval'
contains
  !> \brief Construct an independent physical balance with a prescribed water-temperature root.
  subroutine check_known_root(expected, upstream)
    real(dp), intent(in) :: expected !< Prescribed root [K].
    real(dp), intent(in) :: upstream !< Starting water temperature [K].
    real(dp) :: expected_source, residual, scale
    expected_source = 0.5_dp*area * &
      (river_heat_flux(upstream,ta,sw,lw,albedo,pt,emissivity,kh) + &
       river_heat_flux(expected,ta,sw,lw,albedo,pt,emissivity,kh))
    h0 = water_heat_capacity_volume*q*expected - c2*expected_source
    call solve_river_temperature(q,h0,d,area,upstream,ta,sw,lw,albedo,pt,emissivity,kh,t,source,status,fallback)
    if (status /= 0_i4 .or. abs(t-expected) > 1.0e-6_dp) error stop 'prescribed physical root'
    residual = water_heat_capacity_volume*q*t - h0 - c2*source
    scale = max(abs(water_heat_capacity_volume*q*t),abs(h0),abs(c2*source),1.0_dp)
    if (abs(residual) > 1.0e-10_dp*scale) error stop 'original physical residual'
  end subroutine check_known_root
end program test_river_temperature
