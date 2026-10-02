!> \file    mo_river_temperature.f90
!> \copydoc mo_river_temperature

!> \brief   Open-water river heat balance.
!> \details Solve the open-water heat balance for thermal power relative to 0 K, in watts.
!> \version 0.1
!> \changelog
!! - 2026: v6 open-water scalar heat-balance solver.
!> \authors Sebastian Mueller
!> \date    October 2026
!> \copyright Copyright 2005-\today, the CHS Developers, Sabine Attinger: All rights reserved.
!! This code is released under the LGPLv3+ license \license_note
module mo_river_temperature
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mo_kind, only: dp, i4
  use mo_constants, only: T0_dp, sigma_dp, Psychro_dp, cp_w_dp
  use mo_pet, only: slope_satpressure
  implicit none
  private
  ! cp_w_dp in FORCES is numerically 4.19 kJ kg-1 K-1 despite its J-based comment.
  ! 1000 kg m-3 * 4.19 kJ kg-1 K-1 * 1000 J kJ-1.
  real(dp), parameter, public :: water_heat_capacity_volume = 1000.0_dp * cp_w_dp * 1000.0_dp ! J m-3 K-1
  real(dp), parameter :: min_temperature = 200.0_dp, max_temperature = 400.0_dp
  public :: river_heat_flux, solve_river_temperature

contains
  !> \brief Evaluate atmospheric heat flux into open water in W m-2.
  real(dp) function river_heat_flux(t, ta, shortwave, longwave, albedo, pt, emissivity, kh) result(flux)
    real(dp), intent(in) :: t !< Water temperature [K].
    real(dp), intent(in) :: ta !< Air temperature [K].
    real(dp), intent(in) :: shortwave !< Incoming shortwave [W m-2].
    real(dp), intent(in) :: longwave !< Incoming longwave [W m-2].
    real(dp), intent(in) :: albedo !< Open-water shortwave albedo.
    real(dp), intent(in) :: pt !< Priestley-Taylor coefficient.
    real(dp), intent(in) :: emissivity !< Open-water emissivity.
    real(dp), intent(in) :: kh !< Sensible heat exchange coefficient.
    real(dp) :: r0, delta, eta
    delta = slope_satpressure(ta - T0_dp)
    eta = pt * delta / (Psychro_dp + delta)
    r0 = (1.0_dp - albedo) * shortwave + longwave
    flux = heat_flux_precomputed(t, ta, r0, eta, emissivity, kh)
  end function river_heat_flux

  !> \brief Evaluate heat flux with meteorology-dependent terms computed once for a reach.
  pure real(dp) function heat_flux_precomputed(t, ta, r0, eta, emissivity, kh) result(flux)
    real(dp), intent(in) :: t !< Water temperature [K].
    real(dp), intent(in) :: ta !< Air temperature [K].
    real(dp), intent(in) :: r0 !< Incoming radiation after albedo [W m-2].
    real(dp), intent(in) :: eta !< Priestley-Taylor multiplier.
    real(dp), intent(in) :: emissivity !< Open-water emissivity.
    real(dp), intent(in) :: kh !< Sensible heat exchange coefficient.
    real(dp) :: netrad
    netrad = r0 - emissivity * sigma_dp * t**4
    flux = netrad - eta * max(netrad, 0.0_dp) - kh * (t - ta)
  end function heat_flux_precomputed

  !> \brief Solve the branch-aware heat balance for downstream water temperature.
  !> \details Q is the link outflow in m3 s-1; H0 is the known routed power in W;
  !! D is nu2 times half the surface area in m2. Status 0 is success, 1 no root,
  !! and 2 multiple admissible roots. Scalar root updates are safeguarded by a bracket.
  subroutine solve_river_temperature(q, h0, d, area, tup, ta, shortwave, longwave, albedo, pt, emissivity, kh, &
                                     temperature, surface_source, status, used_fallback)
    real(dp), intent(in) :: q !< Routed link discharge [m3 s-1].
    real(dp), intent(in) :: h0 !< Known Muskingum thermal term [W].
    real(dp), intent(in) :: d !< nu2 times half the physical reach area [m2].
    real(dp), intent(in) :: area !< Physical reach surface area [m2].
    real(dp), intent(in) :: tup !< Incoming water temperature [K].
    real(dp), intent(in) :: ta !< Mean endpoint air temperature [K].
    real(dp), intent(in) :: shortwave !< Mean endpoint shortwave [W m-2].
    real(dp), intent(in) :: longwave !< Mean endpoint longwave [W m-2].
    real(dp), intent(in) :: albedo !< Open-water shortwave albedo.
    real(dp), intent(in) :: pt !< Priestley-Taylor coefficient.
    real(dp), intent(in) :: emissivity !< Open-water emissivity.
    real(dp), intent(in) :: kh !< Sensible heat exchange coefficient.
    real(dp), intent(out) :: temperature !< Solved outgoing water temperature [K].
    real(dp), intent(out) :: surface_source !< Atmospheric source added to link input [W].
    integer(i4), intent(out) :: status !< Zero for one valid root.
    logical, optional, intent(out) :: used_fallback !< Whether bracketed solving was used.
    real(dp) :: r0, eta, delta, tswitch
    if (q <= 0.0_dp .or. area < 0.0_dp) then
      temperature = 0.0_dp
      surface_source = 0.0_dp
      status = 1_i4
      if (present(used_fallback)) used_fallback = .false.
      return
    end if
    delta = slope_satpressure(ta - T0_dp)
    eta = pt * delta / (Psychro_dp + delta)
    r0 = (1.0_dp - albedo) * shortwave + longwave
    tswitch = min_temperature
    if (r0 > 0.0_dp) tswitch = sqrt(sqrt(r0 / (emissivity * sigma_dp)))
    call solve_heat_balance(q, h0, d, area, tup, ta, r0, eta, tswitch, emissivity, kh, &
      temperature, surface_source, status, used_fallback)
  end subroutine solve_river_temperature

  !> \brief Solve the physical balance with forcing coefficients prepared by the caller.
  subroutine solve_heat_balance(q, h0, d, area, tup, ta, r0, eta, tswitch, emissivity, kh, &
                                temperature, surface_source, status, used_fallback)
    real(dp), intent(in) :: q !< Routed discharge [m3 s-1].
    real(dp), intent(in) :: h0 !< Known routed power [W].
    real(dp), intent(in) :: d !< Muskingum coefficient times half the reach area [m2].
    real(dp), intent(in) :: area !< Physical reach surface area [m2].
    real(dp), intent(in) :: tup !< Incoming water temperature [K].
    real(dp), intent(in) :: ta !< Mean endpoint air temperature [K].
    real(dp), intent(in) :: r0 !< Incoming radiation after albedo [W m-2].
    real(dp), intent(in) :: eta !< Priestley-Taylor multiplier.
    real(dp), intent(in) :: tswitch !< Radiation branch boundary [K].
    real(dp), intent(in) :: emissivity !< Open-water emissivity.
    real(dp), intent(in) :: kh !< Sensible heat coefficient [W m-2 K-1].
    real(dp), intent(out) :: temperature !< Solved water temperature [K].
    real(dp), intent(out) :: surface_source !< Atmospheric source [W].
    integer(i4), intent(out) :: status !< Zero for one valid root, one for none, two for multiple.
    logical, optional, intent(out) :: used_fallback !< Whether any bisection was needed.
    real(dp) :: lo, hi, a, a4, a1, a0, b, fup, cvq
    real(dp) :: cuts(3), gx, gl, gu, mid, root, tol, tol_upper, stationary
    real(dp) :: accepted(4)
    integer(i4) :: branch, nc, nroot, j, k
    logical :: valid, bisected

    ! Default to failure and zero outputs so every early return has defined results.
    temperature = 0.0_dp
    surface_source = 0.0_dp
    status = 1_i4
    if (present(used_fallback)) used_fallback = .false.
    if (q <= 0.0_dp .or. area < 0.0_dp) return

    ! These terms are shared by both radiation branches.  The known upstream
    ! surface flux contributes to the constant part of the implicit heat balance.
    cvq = water_heat_capacity_volume * q
    fup = heat_flux_precomputed(tup, ta, r0, eta, emissivity, kh)
    b = h0 + d * fup
    accepted = 0.0_dp
    nroot = 0_i4

    ! The max(netrad,0) latent-heat term gives two quartics: one while
    ! net radiation is positive, and one after it reaches zero.
    do branch = 1_i4, 2_i4
      if (branch == 1_i4) then
        lo = min_temperature
        hi = min(max_temperature, tswitch)
        a = 1.0_dp - eta
      else
        lo = max(min_temperature, tswitch)
        hi = max_temperature
        a = 1.0_dp
      end if
      ! The radiative coefficient is 1-eta while net radiation drives
      ! evaporation; it is 1 when negative radiation switches latent loss off.
      if (hi < lo) cycle

      ! Expand g(T) = Cv*Q*T - H0 - D*(F(Tup)+F(T)) into
      ! a4*T**4 + a1*T + a0 on this branch; only the radiation
      ! coefficient a differs between the sunlit and dark branches.
      a4 = d * a * emissivity * sigma_dp
      a1 = cvq + d * kh
      a0 = -(b + d * (a * r0 + kh * ta))
      ! A quartic can turn once on the branch. Split at its stationary
      ! temperature, when present, so each interval is monotonic and can
      ! be checked for a root by its endpoint residuals.
      cuts = [lo, hi, hi]
      nc = 2_i4
      ! A stationary point solves g'(T) = 4*a4*T**3 + a1 = 0.
      if (abs(a4) > tiny(a4)) then
        stationary = -a1 / (4.0_dp * a4)
        if (stationary > 0.0_dp) then
          stationary = stationary**(1.0_dp / 3.0_dp)
          if (stationary > lo .and. stationary < hi) then
            cuts = [lo, stationary, hi]
            nc = 3_i4
          end if
        end if
      end if
      ! Check each monotonic interval independently. Record roots at cuts
      ! (including a tangent root at a stationary point) before testing for
      ! a sign change in the open interval.
      do j = 1_i4, nc - 1_i4
        call evaluate_balance(cuts(j), gl, tol)
        call evaluate_balance(cuts(j+1), gu, tol_upper)
        if (abs(gl) <= tol) call add_root(cuts(j))
        if (abs(gu) <= tol_upper) call add_root(cuts(j+1))
        if ((gl >= 0.0_dp .and. gu >= 0.0_dp) .or. (gl <= 0.0_dp .and. gu <= 0.0_dp)) cycle
        valid = .false.
        ! Avoid iteration when the quartic term is negligible at every
        ! permitted temperature; solve the remaining linear balance directly.
        if (abs(a4 * max_temperature**4) <= 1.0e-12_dp * &
            max(abs(a1 * max_temperature), abs(a0), 1.0_dp)) then
          ! The quartic contribution is below the residual tolerance: solve the linear limit directly.
          if (abs(a1) > tiny(a1)) then
            root = -a0 / a1
            valid = ieee_is_finite(root)
          end if
        else
          call sparse_quartic_root(a4, a1, a0, tup, cuts(j), cuts(j+1), root, valid, bisected)
          if (present(used_fallback)) used_fallback = used_fallback .or. bisected
        end if
        ! An iterative candidate is provisional: require it to stay on the
        ! interval that produced the bracket and satisfy the original W balance.
        if (valid) valid = root >= cuts(j) .and. root <= cuts(j+1)
        if (valid) then
          call evaluate_balance(root, gx, tol)
          valid = abs(gx) <= tol
        end if
        if (.not.valid) then
          ! Keep the sign-changing interval if the fast solve fails validation.
          ! Bisection stops early when both residual and width tolerances pass;
          ! the final candidate is checked against the physical residual below.
          if (present(used_fallback)) used_fallback = .true.
          lo = cuts(j)
          hi = cuts(j+1)
          do k = 1_i4, 100_i4
            mid = 0.5_dp * (lo + hi)
            call evaluate_balance(mid, gx, tol)
            if (abs(gx) <= tol .and. hi - lo <= 1.0e-6_dp) exit
            if ((gl <= 0.0_dp .and. gx >= 0.0_dp) .or. (gl >= 0.0_dp .and. gx <= 0.0_dp)) then
              hi = mid
            else
              lo = mid
              gl = gx
            end if
          end do
          root = mid
        end if
        call evaluate_balance(root, gx, tol)
        if (abs(gx) <= tol) call add_root(root)
      end do
    end do
    ! Publish only an unambiguous physical root. The atmospheric source is
    ! the trapezoidal mean of the endpoint heat fluxes times full reach area.
    if (nroot == 1_i4) then
      status = 0_i4
      temperature = accepted(1)
      surface_source = area * 0.5_dp * (fup + &
        heat_flux_precomputed(temperature, ta, r0, eta, emissivity, kh))
    else if (nroot > 1_i4) then
      status = 2_i4
    end if
  contains
    !> \brief Evaluate the original heat balance and its relative energy tolerance together.
    subroutine evaluate_balance(t, value, tolerance)
      real(dp), intent(in) :: t !< Trial temperature [K].
      real(dp), intent(out) :: value !< Physical heat-balance residual [W].
      real(dp), intent(out) :: tolerance !< Allowed energy residual [W].
      real(dp) :: power, source
      ! This residual is in watts and matches the routed physical balance;
      ! scale its tolerance to the magnitudes of the competing energy terms.
      power = cvq * t
      source = d * (fup + heat_flux_precomputed(t, ta, r0, eta, emissivity, kh))
      value = power - h0 - source
      tolerance = 1.0e-10_dp * max(abs(power), abs(h0), abs(source), 1.0_dp)
    end subroutine evaluate_balance

    !> \brief Add a unique, physically valid temperature root to the accepted candidates.
    subroutine add_root(t)
      real(dp), intent(in) :: t !< Trial water temperature [K].
      integer(i4) :: m
      if (.not.ieee_is_finite(t)) return
      if (t < min_temperature .or. t > max_temperature) return
      ! Adjacent intervals share cut temperatures; do not count the same root twice.
      do m = 1_i4, nroot
        if (abs(t - accepted(m)) < 1.0e-6_dp) return
      end do
      if (nroot < size(accepted)) then
        nroot = nroot + 1_i4
        accepted(nroot) = t
      end if
    end subroutine add_root
  end subroutine solve_heat_balance

  !> \brief Find a root of a4*T**4+a1*T+a0 on one monotonic branch without allocating arrays.
  subroutine sparse_quartic_root(a4, a1, a0, seed, lower, upper, root, valid, bisected)
    real(dp), intent(in) :: a4 !< Quartic coefficient.
    real(dp), intent(in) :: a1 !< Linear coefficient.
    real(dp), intent(in) :: a0 !< Constant coefficient.
    real(dp), intent(in) :: seed !< Incoming water-temperature seed [K].
    real(dp), intent(in) :: lower !< Lower branch endpoint [K].
    real(dp), intent(in) :: upper !< Upper branch endpoint [K].
    real(dp), intent(out) :: root !< Candidate root [K].
    logical, intent(out) :: valid !< Whether the normalized polynomial converged.
    logical, intent(out) :: bisected !< Whether a safeguarded midpoint was needed.
    real(dp) :: c4, c1, c0, scale, x, x2, x3, f, df, step, candidate, lo, hi, fl, checkpoint
    integer(i4) :: iteration

    ! Solve in x = T/300 K so quartic and linear terms have comparable scales.
    c4 = a4 * 300.0_dp**4
    c1 = a1 * 300.0_dp
    scale = max(abs(c4), abs(c1), abs(a0))
    ! Keep the seed as the failure result; the caller can reject it via valid.
    root = seed
    valid = .false.
    bisected = .false.
    if (.not.ieee_is_finite(scale) .or. scale <= 0.0_dp) return
    c4 = c4 / scale
    c1 = c1 / scale
    c0 = a0 / scale
    lo = lower / 300.0_dp
    hi = upper / 300.0_dp
    ! Clamp the upstream-temperature seed to this monotonic branch's bracket.
    x = min(max(seed / 300.0_dp, lo), hi)
    ! The branch endpoints bracket at most one root; fl tracks the lower endpoint.
    fl = (c4 * lo**3 + c1) * lo + c0
    checkpoint = hi - lo
    do iteration = 1_i4, 50_i4
      ! Reuse powers for both polynomial value and its Newton derivative.
      x2 = x * x
      x3 = x2 * x
      f = (c4 * x3 + c1) * x + c0
      df = 4.0_dp * c4 * x3 + c1
      if (abs(f) < 1.0e-12_dp) then
        root = x * 300.0_dp
        valid = .true.
        return
      end if
      if ((fl <= 0.0_dp .and. f <= 0.0_dp) .or. (fl >= 0.0_dp .and. f >= 0.0_dp)) then
        lo = x
        fl = f
      else
        hi = x
      end if
      ! Newton is fast near a root; the bracket midpoint remains the safe default.
      candidate = 0.5_dp * (lo + hi)
      if (abs(df) > epsilon(df) * max(abs(4.0_dp*c4*x3), abs(c1), tiny(df))) then
        step = f / df
        candidate = x - step
      else
        bisected = .true.
      end if
      ! Periodically enforce bracket progress even if all proposed root updates stay inside it.
      if (.not.ieee_is_finite(candidate) .or. candidate <= lo .or. candidate >= hi .or. candidate == x) then
        candidate = 0.5_dp * (lo + hi)
        bisected = .true.
      else if (mod(iteration, 4_i4) == 0_i4) then
        if (hi - lo > 0.5_dp * checkpoint) then
          candidate = 0.5_dp * (lo + hi)
          bisected = .true.
        end if
        checkpoint = hi - lo
      end if
      x = candidate
    end do
    ! Return the last estimate for diagnostics, but leave valid false on exhaustion.
    root = x * 300.0_dp
  end subroutine sparse_quartic_root
end module mo_river_temperature
