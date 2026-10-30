!> \file mo_lake_water_balance.F90
!> \copydoc mo_lake_water_balance

!> \brief simple lake water balance module

!> \details  lake water balance module for balancing lake inflow from mRM,
!>            precipitation, evaporation, percolation and spill. The solution
!>            gives the lake states i.e. volume, area and level for each time step.

!> \authors Pallav K Shrestha, Luis Samaniego, Rohini Kumar

!> \date March 2018

!> \version 0.1
!> \changelog
!! - October 2026: v6 component ownership and restart integration; retained donor arithmetic.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mLM
MODULE mo_lake_water_balance

  USE mo_kind,   only : i4, i8, dp
  USE mo_mlm_constants, only : error_lake_vol, nIterations_lake
  IMPLICIT NONE
  private


  PUBLIC :: lake_balance     ! calculates lake water balance to give lake states (V, A, h)
  PUBLIC :: interpolate_lake_states

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         lake_balance

  !     PURPOSE
  !>        \brief Perform lake water balance

  !>        \details Perform lake water balance between inlets, storage and outlet at each time step

  !     CALLING SEQUENCE

  !     RETURN
  !         None

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Mar 2018

  !     MODIFICATIONS
  !        Pallav Shrestha Sep 2018 - Moved the volume area elevation (EVA) read procedure to
  !                               mo_mlm_read_config as a sub routine "read_lake_bathymetry"
  !        Pallav Shrestha Nov 2019 - Added procedure for adjusting volume when reservoir overtopps/ gets empty
  !        Pallav Shrestha Jun 2022 - Added IBT
  !        Pallav Shrestha Oct 2026 - mLM v6 migration: ported from mo_mlm_balance.
  !                               (1) SI units: the v5 water-balance arithmetic carried area in
  !                                   km2 and volume in mcm (the +-1e6 factors); these are
  !                                   dropped so the balance runs in m2 / m3 throughout. The
  !                                   closure tolerance is correspondingly 1 m3 (see module head).
  !                               (2) Nomenclature: h/A/V -> lake_level/lake_surface_area/
  !                                   lake_volume; the hAV curve is passed as three explicit
  !                                   columns (hav_lake_level/_surface_area/_volume) instead of
  !                                   lake_bathy(:,1:3); DCL/DCL_A/DCL_V ->
  !                                   lake_spill_reference_level / lake_spill_reference_surface_area_bathymetry_derived
  !                                   / lake_spill_reference_volume_bathymetry_derived.
  !                               (3) State time levels: the v5 2-element (prev,cur) arrays are
  !                                   replaced by explicit *_previous (in) and current (out)
  !                                   scalar arguments, following the v6 mRM convention.
  !                               (4) No globals: tsHydroSecs is now an argument (was
  !                                   DaySecs / nTStepDay); the unused mo_utils / mo_constants /
  !                                   mo_common_mHM_mRM_variables imports and the unused tt
  !                                   argument were dropped. Physics/equations are unchanged.
  !        Pallav Shrestha Oct 2026 - mLM v6 migration (Phase 5): the solver constants
  !                                   error_lake_vol and nIterations_lake were moved to
  !                                   mo_mlm_constants (now that it exists) and are used from
  !                                   there; values unchanged.



  recursive SUBROUTINE lake_balance(naturalFlag, pr, evap, inflow, perc, wateruse, spill, ibt, &
                          hav_lake_level, hav_lake_surface_area, hav_lake_volume, &
                          lake_spill_reference_level, lake_spill_reference_surface_area_bathymetry_derived, &
                          lake_spill_reference_volume_bathymetry_derived, &
                          lake_level_previous, lake_surface_area_previous, lake_volume_previous, &
                          lake_level, lake_surface_area, lake_volume, volume_adjustment, tsHydroSecs)

    IMPLICIT NONE

    logical,    intent(in) :: naturalFlag !< [-] regulation flag (regulated vs natural lake)

    real(dp), intent(IN) :: pr !< [mm.TS-1]  disaggregated lake precipitation
    real(dp), intent(IN) :: evap !< [mm.TS-1]  disaggregated lake evaporation
    real(dp), intent(IN) :: inflow !< [m3.s-1]   lake inflow from mRM
    real(dp), intent(IN) :: perc !< [mm.TS-1]  lake percolation [per routing time step]
    real(dp), intent(IN) :: wateruse !< [m3.s-1]   lake wateruse
    real(dp), intent(IN) :: spill !< [m3.s-1]   lake spill
    real(dp), intent(IN) :: ibt !< [m3.s-1]   lake inter-basin transfer (IBT)

    real(dp), intent(in) :: hav_lake_level(:) !< hAV level column [m]
    real(dp), intent(in) :: hav_lake_surface_area(:) !< hAV area column [m2]
    real(dp), intent(in) :: hav_lake_volume(:) !< hAV volume column [m3]
    real(dp), intent(IN) :: lake_spill_reference_level !< [masl] dam crest level/ outflow elevation for natural lakes
    real(dp), intent(IN) :: lake_spill_reference_surface_area_bathymetry_derived !< [m2]   area corresponding to the crest
    real(dp), intent(IN) :: lake_spill_reference_volume_bathymetry_derived !< [m3]   volume corresponding to the crest

    real(dp), intent(IN) :: lake_level_previous !< [masl]  previous-step lake level
    real(dp), intent(IN) :: lake_surface_area_previous !< [m2]    previous-step lake surface area
    real(dp), intent(IN) :: lake_volume_previous !< [m3]    previous-step lake volume

    real(dp), intent(OUT) :: lake_level !< [masl]  current lake level
    real(dp), intent(OUT) :: lake_surface_area !< [m2]    current lake surface area
    real(dp), intent(OUT) :: lake_volume !< [m3]    current lake volume

    real(dp), intent(OUT) :: volume_adjustment !< [m3]    volume adjusted for the current time step

    real(dp), intent(IN) :: tsHydroSecs !< [s]     number of seconds in a hydrological time step

    !! local vars

    !--------------------------------------------------------------------------------------
    ! Usage structure of local variables in the water balance (WB) iteration
    !--------------------------------------------------------------------------------------
    !
    !     delta_V = (dv =) V(t) - V(t-1)
    !
    !                                     where t is the current time step
    !
    !       Steps -
    !         1. A delta_V_trial is assumed to get V and then A and h from EVA curves
    !            for the current time step
    !         2. Based on this trial set of V, A, h, delta_V is obtained from WB equation
    !
    !     delta_V = [pr(t) - evap(t)] * [A(t-1) + A(t)]/2 - perc(t) * A(t-1) + [inflow(t)
    !                 - wateruse(t) - spill(t) + ibt(t) ] * tsHydroSecs
    !
    !            Note: seepage perc(t) is removed over the previous-step area A(t-1) only, so the lake
    !                  loss equals the volume added to the groundwater store (water conservation);
    !                  the P - E term uses the average area as before.
    !
    !            Note: wateruse(t) = sum[hydelflow(t), enviflow(t), watersupply(t), irrigation(t)]
    !
    !                  tsHydroSecs is number of seconds in a hydrological time step. This converts flows to
    !                  volume for the hydrological time step.
    !
    !         3. The difference between the two values of dV is the volume error.
    !
    !     volume error = dV from WB equation - dV assumed
    !
    !         4. If ...
    !                   volume error is less than allowed error then the iteration stops
    !                   and the lake states for the current time step is finalized.
    !            Else...
    !                   the next delta_V_trial is calculated from residuals of previous
    !                   two iterations using "SECANT METHOD" and the next iteration is
    !                   carried out
    !
    !                                                  dV_trial(n-1) - dV_trial(n-2)
    !     dV_trial(n) = dV_trial(n-1) - Verror(n-1) x  -----------------------------
    !                                                    Verror(n-1) - Verror(n-2)
    !
    !-----------------------------------------------------------------------------------------

    real(dp)       :: delta_V_trial_n, &   ! dV assumed for current (n)th iteration
                      delta_V_trial_n_1, & ! dV assumed for (n-1)th iteration
                      delta_V_trial_n_2, & ! dV assumed for (n-2)th iteration

                      !
                      delta_V_n, &         ! dV from WB equation for current (n)th iteration
                      delta_V_n_1, &       ! dV from WB equation for (n-1)th iteration
                      delta_V_n_2, &       ! dV from WB equation for (n-2)th iteration


                      Verror_n, &          ! volume error for current (n)th iteration
                      Verror_n_1, &        ! volume error for (n-1)th iteration
                      Verror_n_2, &        ! volume error for (n-2)th iteration

                      h_n_1, &             ! water level for (n-1)th iteration
                      h_n_2, &             ! water level for (n-2)th iteration
                      A_n_1, &             ! surface area for (n-1)th iteration
                      A_n_2, &             ! surface area for (n-2)th iteration
                      V_n_1, &             ! volume for (n-1)th iteration
                      V_n_2                ! volume for (n-2)th iteration

    integer(i4)    :: iIteration


      !---------------------------------------------------------
      ! calculating lake states
      !---------------------------------------------------------

      ! Calculating the first two iteration error functions
      ! Assumption: The water level fluctuates by -2cm/+1cm

      h_n_1 = lake_level_previous - 0.02_dp
      h_n_2 = lake_level_previous + 0.01_dp

      call interpolate_lake_states(.not. naturalFlag, hav_lake_level, hav_lake_volume, &
        lake_spill_reference_level, lake_spill_reference_volume_bathymetry_derived, h_n_1, V_n_1)
      call interpolate_lake_states(.not. naturalFlag, hav_lake_level, hav_lake_volume, &
        lake_spill_reference_level, lake_spill_reference_volume_bathymetry_derived, h_n_2, V_n_2)

      delta_V_trial_n_2 = V_n_2 - lake_volume_previous
      delta_V_trial_n_1 = lake_volume_previous - V_n_1

      call interpolate_lake_states(.not. naturalFlag, hav_lake_level, hav_lake_surface_area, &
        lake_spill_reference_level, lake_spill_reference_surface_area_bathymetry_derived, h_n_2, A_n_2)
      call interpolate_lake_states(.not. naturalFlag, hav_lake_level, hav_lake_surface_area, &
        lake_spill_reference_level, lake_spill_reference_surface_area_bathymetry_derived, h_n_1, A_n_1)


      delta_V_n_2 = ( pr - evap ) / 1000.0_dp * 0.5_dp * ( lake_surface_area_previous + A_n_2 ) &
       - perc / 1000.0_dp * lake_surface_area_previous &
       + ( inflow - wateruse - spill + ibt ) * tsHydroSecs
      delta_V_n_1 = ( pr - evap ) / 1000.0_dp * 0.5_dp * ( lake_surface_area_previous + A_n_1 ) &
       - perc / 1000.0_dp * lake_surface_area_previous &
       + ( inflow - wateruse - spill + ibt ) * tsHydroSecs

          ! where -
                  ! Met observations (mm) divided by 1000 converts them to m.hydrologicalTS-1 ;
                  ! multiplied by the average surface area (m2) gives the atmospheric volume (m3).

                  ! Flow rates (m3.s-1) multiplied by number of seconds per hydrological time step
                  ! (tsHydroSecs) changes them to m3.hydroTS-1 .

                  ! The whole expression gives delta_V directly in m3 (SI). In v5 the area was in
                  ! km2 and the volume in mcm, hence the +-1e6 factors that are dropped here.

      Verror_n_2 = delta_V_n_2 - delta_V_trial_n_2
      Verror_n_1 = delta_V_n_1 - delta_V_trial_n_1

      volume_adjustment = 0._dp ! volume adjusted for avoiding reservoir volume instability issues in this time step

      ! Calculation for third iteration onwards
      IterationLoop: do iIteration = 1, nIterations_lake

        ! SECANT METHOD formula for next iteration volume difference estimate
        delta_V_trial_n = delta_V_trial_n_1 - Verror_n_1 * (( delta_V_trial_n_1 &
          - delta_V_trial_n_2) / ( Verror_n_1 - Verror_n_2))

        lake_volume = lake_volume_previous + delta_V_trial_n

        ! Is the RESERVOIR EMPTY??
        if (lake_volume < 0) then

          ! Reservoir is empty
          lake_surface_area = 0._dp
          lake_level = hav_lake_level(1)

          delta_V_n = ( pr - evap ) / 1000.0_dp * 0.5_dp * (lake_surface_area_previous + lake_surface_area ) &
            - perc / 1000.0_dp * lake_surface_area_previous &
            + ( inflow - wateruse - spill + ibt ) * tsHydroSecs

          lake_volume = lake_volume_previous + delta_V_n ! This is considered as the final value, and is negative

          volume_adjustment = -lake_volume ! To store +ve value. +ve: Reservoir empty, -ve: dam overtop,

          lake_volume = 0._dp

        else

          call interpolate_lake_states(.not. naturalFlag, hav_lake_volume, hav_lake_level, &
            lake_spill_reference_volume_bathymetry_derived, lake_spill_reference_level, lake_volume, lake_level, volume_adjustment)
                                                       ! volume_adjustment gives the volume
                                                       ! adjusted to prevent reservoir volume instabilities
          call interpolate_lake_states(.not. naturalFlag, hav_lake_level, hav_lake_surface_area, &
            lake_spill_reference_level, lake_spill_reference_surface_area_bathymetry_derived, lake_level, lake_surface_area)

          delta_V_n = ( pr - evap ) / 1000.0_dp * 0.5_dp * (lake_surface_area_previous + lake_surface_area ) &
            - perc / 1000.0_dp * lake_surface_area_previous &
            + ( inflow - wateruse - spill + ibt ) * tsHydroSecs

          ! Percolation is intentionally applied over the previous-step lake area (A_prev), not the
          ! average area used for the net atmospheric flux (P - E). perc is computed once (as a depth)
          ! before the water balance, and the same perc/1000 * A_prev volume is removed here and added
          ! to the groundwater store (lake_gws_baseflow), so the lake->groundwater transfer conserves
          ! water exactly. (The P - E term keeps the average area because precipitation and evaporation
          ! act on the dynamic water surface over the step.) - PKS Oct 2018; conservation closed 2026

          Verror_n = delta_V_n - delta_V_trial_n

        end if

        ! if the volume has been adjusted for reservoir volume instabilities, the corresponding
        ! h, A, V are the states for this iteration. Further iterations halted!
        if (volume_adjustment /= 0._dp) exit


        if ( abs(Verror_n) .le. error_lake_vol) then
          exit
        else
          delta_V_trial_n_2 = delta_V_trial_n_1
          delta_V_trial_n_1 = delta_V_trial_n
          Verror_n_2   =  Verror_n_1
          Verror_n_1   =  Verror_n
        endif

      end do IterationLoop

  END SUBROUTINE lake_balance


  ! ------------------------------------------------------------------
  !     NAME
  !         interpolate_lake_states

  !     PURPOSE
  !>        \brief linear interpolation of lake states

  !>        \details linear interpolation of lake states based on EVA curves

  !     CALLING SEQUENCE

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         Both series used for interpolation should span the whole interpolation length
  !         (the v6 hAV curves are exactly sized, so no nodata padding is present).

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Mar 2018

  !     MODIFICATIONS
  !        Pallav Shrestha Oct 2026 - mLM v6 migration: dropped the unused optimize import and
  !                               tstep argument; the source-series length is now size(series1)
  !                               (the v6 hAV curves carry no nodata padding, unlike the v5
  !                               fixed-length bathy arrays), removing the mo_common_constants
  !                               nodata_dp dependency. Interpolation logic is unchanged.

  recursive SUBROUTINE interpolate_lake_states(flag, series1, series2, max1, max2, val1, val2, vol_adj)

    IMPLICIT NONE

    logical,  intent(in) :: flag !< Flag for anthropogenic regulation i.e. reservoir (TRUE), natural lake (FALSE)
    real(dp), dimension(:), intent(in) :: series1 !< The first/ source array
    real(dp), dimension(:), intent(in) :: series2 !< The second/ destination array
    real(dp), intent(in) :: max1 !< crest (DCL) value for the source lake state
    real(dp), intent(in) :: max2 !< crest (DCL) value for the target lake state
    real(dp), intent(inout) :: val1 !< Source value
    real(dp), intent(out) :: val2 !< Value to be interpolated
    real(dp), optional, intent(out) :: vol_adj !< Volume adjustment for reservoir instability conditions.
                                                      ! When this argument is present, source series is volume

    ! local vars
    real(dp)              :: val1_minus               ! Lower neighbor of val1
    real(dp)              :: val2_minus               ! Corresponding value of val1_minus in second array
    real(dp)              :: val1_plus                ! Higher neighbor of val1
    real(dp)              :: val2_plus                ! Corresponding value of val1_plus in second array
    logical               :: interpolate
    integer(i4)           :: ssize1                   ! size of the input series
    integer(i4), dimension(2) :: pos                  ! position index in source array (series1)

    if (present(vol_adj)) vol_adj = 0._dp
    ssize1 = size(series1)


    if ( val1 <= series1(1) ) then

      ! hook interpolated value to the lowest value
      val1 = series1(1)
      val2 = series2(1)

    else

      interpolate = .FALSE.

      ! get the position of series 1 value which is lower/ equal to val1
      pos(1:1) = maxloc(series1(1:ssize1), series1(1:ssize1) <= val1)

      ! reservoir or natural lake?
      if (flag) then
        Reservoir: if ( val1 <= max1 ) then ! interpolate
          interpolate = .TRUE.
        else ! no interpolation beyond DCL values
          ! Beyond DCL, Dam --  OVERTOPPED! Spill excess water and bring water level to DCL
          ! If needed, estimate volume adjusted i.e. water removed for spilling
          if (present(vol_adj)) vol_adj = max1 - val1
                                          ! When this argument is present, source series is volume
                                          ! Negative value represents "removal" of water
          ! hook source and interpolated values to their corresponding maximum
          val1 = max1
          val2 = max2
        end if Reservoir
      else
        ! Natural lake
        interpolate = .TRUE.
      endif


      if ( interpolate ) then

        if ( val1 >= series1(ssize1) ) then
          ! extrapolate
          pos(2) = pos(1)
          pos(1) = pos(1) - 1
        else
          ! interpolate
          pos(2) = pos(1) + 1
        end if

        ! get the neighboring values for the two series
        val1_minus  = series1( pos(1) )
        val1_plus   = series1( pos(2) )
        val2_minus  = series2( pos(1) )
        val2_plus   = series2( pos(2) )

        ! interpolate
        val2 = val2_minus + ( val1 - val1_minus) * (  val2_plus - val2_minus)/(val1_plus - val1_minus)

      end if

    end if


  END SUBROUTINE interpolate_lake_states

  !!---------------------------------------------------------
END MODULE mo_lake_water_balance
