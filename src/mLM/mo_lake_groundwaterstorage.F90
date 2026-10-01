!> \file mo_lake_groundwaterstorage.f90

!> \brief groundwater storage module

!> \details  lake groundwater storage module for calculating the percolation loss &
!>            the baseflow transmission to the river downstream of the lake.
!>            The percolation is based on the lake percolation coefficient
!>            and the head above the bed level. The percolated water adds to
!>            the groundwater storage of the reservoir. Then, the baseflow
!>            at each time step is a fraction which is the product of the
!>            baseflow coefficient and the GW storage

!> \authors Pallav K Shrestha, Luis Samaniego

!> \date October 2018

! Modifications:


MODULE mo_lake_groundwaterstorage

  USE mo_kind,   only : i4, i8, dp
  IMPLICIT NONE


  PUBLIC :: lake_seepage     ! calculates the water seepage/ percolation beneath the reservoir
  PUBLIC :: lake_gws_baseflow! updates GWS and calculates the baseflow back to the river downstream from the reservoir.


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         lake_seepage

  !     PURPOSE
  !>        \brief Calculates lake seepage and baseflow

  !>        \details Calculates the percolation loss &
  !>                 the baseflow transmission to the river downstream of the lake.
  !>                 The percolation is based on the lake percolation coefficient
  !>                 and the head above the bed level. The percolated water adds to
  !>                 the groundwater reservoir of the reservoir. Then, the baseflow
  !>                 at each time step is a fraction which is the product of the
  !>                 baseflow coefficient and the GW storage. The
  !>                 baseflow is deducted from the groundwater storage of the reservoir.

  !     CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "logical     :: naturalFlag"                           regulation flag
  !>       \param[in] "real(dp), dimension(:) :: hav_lake_level"          [masl]  hAV curve level column
  !>       \param[in] "real(dp), dimension(:) :: hav_lake_surface_area"   [m2]    hAV curve surface area column
  !>       \param[in] "real(dp) :: lake_max_level_reported"                    [masl] dam crest level (DCL)
  !>       \param[in] "real(dp) :: lake_max_surface_area_bathymetry_derived"   [m2]   area at the crest (DCL_A)
  !>       \param[in] "real(dp) :: lake_level"                            [masl]  lake level
  !>       \param[in] "real(dp) :: percolation_coefficient"               [TS-1] percolation coefficient for lake (per hydrological time step)
  !>       \param[in] "real(dp) :: exponent_for_percolation"              [-]    exponential constant that controls the
  !                                                                              fraction kfrac apart for estimating the
  !                                                                              available water for percolation

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>       \param[out] "real(dp) :: perc"                         [mm TS-1]  lake percolation

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Oct 2018

  !     MODIFICATIONS
  !        Pallav Shrestha Oct 2026 - mLM v6 migration: ported from mo_mlm_groundwaterstorage.
  !                               (1) SI units: the hAV surface areas are in m2 (were km2); the
  !                                   area ratio wettedArea/upperArea is dimensionless either way,
  !                                   and the head-to-depth factor (* 1000) converts m -> mm, so the
  !                                   percolation stays in mm per hydrological time step. Physics
  !                                   unchanged.
  !                               (2) Nomenclature: h -> lake_level; the hAV curve is passed as the
  !                                   two explicit columns hav_lake_level / hav_lake_surface_area
  !                                   (were lake_bathy(:,1:2)); DCL/DCL_A ->
  !                                   lake_max_level_reported / lake_max_surface_area_bathymetry_derived;
  !                                   kp -> percolation_coefficient; exp_chi -> exponent_for_percolation.
  !                               (3) No globals: interpolate_lake_states is used from
  !                                   mo_lake_water_balance; the source-series length is size()
  !                                   (the v6 hAV curves carry no nodata padding), removing the
  !                                   mo_common_constants nodata_dp dependency; the unused tt
  !                                   argument is dropped. lake_level is intent(IN) (a local working
  !                                   copy is passed to the interpolate call, which takes it inout).
  !                               Note: percolation_coefficient is the per-hydrological-time-step
  !                               value; the container forms it from the per-day lakes_1 parameter
  !                               as (tsHydroSecs/86400) / percolation_coefficient, reproducing the
  !                               v5 c2TSTu/kPerco form.



  SUBROUTINE lake_seepage(naturalFlag, hav_lake_level, hav_lake_surface_area, &
                          lake_max_level_reported, lake_max_surface_area_bathymetry_derived, &
                          lake_level, percolation_coefficient, exponent_for_percolation, perc)

    USE mo_lake_water_balance,        ONLY: interpolate_lake_states


    IMPLICIT NONE

    logical,    intent(in):: naturalFlag         ! [-] regulation flag (regulated vs natural lake)

    real(dp), dimension(:),  intent(IN) &
                          ::  hav_lake_level, &        ! [masl]  hAV curve level column
                              hav_lake_surface_area    ! [m2]    hAV curve surface area column

    real(dp), intent(IN)  :: lake_max_level_reported                   ! [masl]  dam crest level/ outflow elevation for natural lakes
    real(dp), intent(IN)  :: lake_max_surface_area_bathymetry_derived  ! [m2]    area corresponding to the crest

    real(dp),    intent(IN) &
                          ::  lake_level     ! [masl]  lake level

    real(dp),                   intent(IN) &
                          ::  percolation_coefficient, & ! [TS-1]  percolation coefficient for lake (per hydrological time step)
                              exponent_for_percolation   ! [-]     exponential constant that controls the fraction kfrac apart for estimating
                                             !         available water for percolation

    real(dp),                   intent(OUT) &
                          :: perc            ! [mm TS-1]  lake percolation


    !! local vars
    integer(i4)           :: iEntry, nlines
    real(dp)              :: upper_level, &  ! [masl]   working copy of the lake level for the interpolation call
                             avgBedlevel, &  ! [masl]   average bedlevel of lake slice
                             wettedArea, &   ! [m2]     wetted surface area of lake slice
                                             !          assumed to be ~ lake surface area
                                             !          of lake slice
                             upperArea, &    ! [m2]     The current lake surface area
                             kfrac           ! [-]      fractional multiplier to get available
                                             !          water for lake percolation

    !----------------------------------------------------------------------------
    !
    ! ================== LAKE/ RESERVOIR SEEPAGE ================================
    !
    !     Percolation:
    !
    !         p = ∑(kp x kfrac x H x As/A)
    !
    !
    !         where, p : percol, percolation from the lake bed slice [mm TS-1]
    !                kp: percolation coefficient [TS-1]
    !             kfrac: fractional multiplier to get available water for lake percolation [-]
    !
    !                     kfrac = h / e^(a)
    !
    !                H : pressure head at the lake bed slice [mm]
    !                As: area of the lake bed slice [m2]
    !                A : total area of the lake bed (assumed to be the lake surface area) [m2]
    !
    !
    !     \                                 .___
    !      \--------------------------------|   |
    !       \ slice n                       |   |
    !       |\------------------------------|    \
    !       v \ slice .. ___level___________|     \
    !          \----------------------- ^ --|      \
    !     P     \ slice ..              |   |       \
    !      E    |\--------------------- H --|        \
    !       R   v \ slice 3             E   |         \
    !        C     \------------------- A --|  DAM     \
    !         O     \ slice 2           D   |           \
    !          L    |\----------------- | --|            \
    !           A   v \ slice 1         v   |             \
    !            T     \--------------------'--------------'---.
    !             I       |       |      |                      \
    !              O      v       v      v                       \________.   downstream water
    !               N                                 baseflow             \--------
    !                     Groundwater Storage -- -- -- -- -- -- -- -->      \_______
    !
    !------------------------------------------------------------------------

    !---------------------------------
    ! Lake percolation (slice wise)
    !---------------------------------

    ! interpolating the current lake surface area: "upperArea"
    ! (a local working copy of the lake level is passed; interpolate_lake_states takes it inout)
    upper_level = lake_level
    call interpolate_lake_states(.not. naturalFlag, hav_lake_level, hav_lake_surface_area, &
      lake_max_level_reported, lake_max_surface_area_bathymetry_derived, upper_level, upperArea)

    if (upperArea == 0._dp) then
      ! the reservoir is empty, ie zero seepage
      perc = 0._dp
    else

      ! Slice 1 ...
        ! The average bed level of first slice is average of slice's top level
        ! and deepest level of lake. Since deepest level info is missing,
        ! its assumed here that it is equal to the bed level of first slice.

        avgBedlevel = hav_lake_level(1)
        wettedArea = hav_lake_surface_area(1)

        ! All the head will not percolate. We need a fraction that is a function of head
        kfrac = (lake_level - avgBedlevel) / exp(exponent_for_percolation)

                      ! kfrac = fraction for getting water available for percolation.

                      ! Larger the head, larger should be the water available for percolation. The above
                      ! expression is a linear relationship to get the influence of head on the
                      ! water available for percolation.

                      ! exp(exponent_for_percolation) must have minimum value of the head itself [max head = dam crest]
                      ! and can be many folds larger to other end of the range.

                      ! Since exponent_for_percolation controls the fraction kfrac apart from the head, it is a parameter
                      ! related to the permeability of soil at lake bed which could be regionalized based on
                      ! soil later - PKS

                      ! exponent_for_percolation [-] is an exponential constant and a model parameter [external for now - PKS].
                      ! kfrac [-] a function of the head and assumed to be dimensionless.

        ! percolation from bottom most slice (mm TS-1)
        perc = percolation_coefficient * kfrac * (lake_level - avgBedlevel) * 1000.0_dp * (wettedArea / upperArea)
                      ! kfrac is a dimentionless fraction which when multiplied to the head gives
                      ! the available water for lake percolation. This available water is then multiplied
                      ! by the recession coefficient to get percolation.

                      ! The weighted wetted area helps in getting better representative percolation water
                      ! (mm) based on depth and wetted area at each reservoir slice.

        nlines = size(hav_lake_level)

      ! Slice 2 onwards ...
      EntriesLoop: do iEntry = 2, nlines

        if (lake_level < hav_lake_level(iEntry)) then

          ! average height of current entry
          avgBedlevel = 0.5_dp * (lake_level + hav_lake_level(iEntry - 1))

          ! additional area of current entry
          wettedArea = upperArea - hav_lake_surface_area(iEntry - 1)

          exit

        else

          ! average height of current entry
          avgBedlevel = 0.5_dp * (hav_lake_level(iEntry) + hav_lake_level(iEntry - 1))

          ! additional area of current entry
          wettedArea = hav_lake_surface_area(iEntry) - hav_lake_surface_area(iEntry - 1)

        end if

        ! All the head will not percolate. We need a fraction that is a function of head
        kfrac = (lake_level - avgBedlevel) / exp(exponent_for_percolation)

        ! percolation from this slice
        perc = perc + percolation_coefficient * kfrac * (lake_level - avgBedlevel) &
                * 1000.0_dp * (wettedArea / upperArea)          ![mm "hydrological" TS-1]

      end do EntriesLoop

    end if


  END SUBROUTINE lake_seepage

  ! ------------------------------------------------------------------

  !     NAME
  !         lake_gws_baseflow

  !     PURPOSE
  !>        \brief Calculates lake GWS and baseflow

  !>        \details Calculates lake GWS and baseflow

  !     CALLING SEQUENCE

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: perc"                           [mm TS-1]  lake percolation to groundwater storage
  !>       \param[in] "real(dp) :: lake_surface_area"              [m2]       lake surface area
  !>       \param[in] "real(dp) :: lake_gw_storage_previous"       [m3]       previous-step lake groundwater storage
  !>       \param[in] "real(dp) :: baseflow_recession_coefficient" [TS-1]     groundwater storage recession coefficient (per hydrological time step)
  !>       \param[in] "real(dp) :: tsHydroSecs"                    [s]        number of seconds in a hydrological time step

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>       \param[out] "real(dp) :: lake_gw_storage"               [m3]       current lake groundwater storage
  !>       \param[out] "real(dp) :: baseflow"                      [m3 s-1]   baseflow from groundwater storage

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Feb 2022

  !     MODIFICATIONS
  !        Pallav Shrestha Oct 2026 - mLM v6 migration: ported from mo_mlm_groundwaterstorage.
  !                               (1) SI units: the groundwater storage is carried in m3 (was mcm)
  !                                   and the percolation volume perc/1000 * lake_surface_area is in
  !                                   m3 directly (lake_surface_area is m2); the final mcm -> m3
  !                                   factor (* 1e6) on the baseflow is therefore dropped, leaving
  !                                   baseflow = baseflow / tsHydroSecs in m3 s-1. Physics unchanged.
  !                               (2) Nomenclature: percol -> perc; A -> lake_surface_area; kbf ->
  !                                   baseflow_recession_coefficient. The v5 intent(IN) comment
  !                                   "percol [mcm]" was a mislabel: the body consumed it as mm
  !                                   (perc/1000 * area), which is what is carried here.
  !                               (3) State time levels: the v5 2-element satSTW(prev,cur) array is
  !                                   replaced by explicit lake_gw_storage_previous (in) and
  !                                   lake_gw_storage (out) scalar arguments, following the v6 mRM
  !                                   convention.
  !                               Note: baseflow_recession_coefficient is the per-hydrological-time-step
  !                               value; the container forms it from the per-day lakes_1 parameter as
  !                               (tsHydroSecs/86400) / baseflow_recession_coefficient, reproducing the
  !                               v5 c2TSTu/kBaseFlow form.



  SUBROUTINE lake_gws_baseflow( perc, lake_surface_area, lake_gw_storage_previous, lake_gw_storage, &
                                baseflow_recession_coefficient, tsHydroSecs, baseflow )


    IMPLICIT NONE

    real(dp),   intent(IN)  :: perc            ! [mm TS-1] lake percolation to groundwater storage
    real(dp),   intent(IN)  :: lake_surface_area  ! [m2]  lake surface area
    real(dp),   intent(IN)  :: lake_gw_storage_previous  ! [m3]  previous-step lake groundwater storage
    real(dp),   intent(OUT) :: lake_gw_storage ! [m3]    current lake groundwater storage
    real(dp),   intent(IN)  :: baseflow_recession_coefficient  ! [TS-1]  groundwater storage recession coefficient (per hydrological time step)
    real(dp),   intent(IN)  :: tsHydroSecs     ! [s]    number of seconds in a hydrological timestep
    real(dp),   intent(OUT) :: baseflow        ! [m3.s-1] baseflow from ground water storage

    !----------------------------------------------------------------------------
    !
    ! ================== LAKE/ RESERVOIR GROUNDWATER STORAGE and BASEFLOW ================================
    !
    !
    !     \                                 .___
    !      \--------------------------------|   |
    !       \ slice n                       |   |
    !       |\------------------------------|    \
    !       v \ slice .. ___level___________|     \
    !          \----------------------- ^ --|      \
    !     P     \ slice ..              |   |       \
    !      E    |\--------------------- H --|        \
    !       R   v \ slice 3             E   |         \
    !        C     \------------------- A --|  DAM     \
    !         O     \ slice 2           D   |           \
    !          L    |\----------------- | --|            \
    !           A   v \ slice 1         v   |             \
    !            T     \--------------------'--------------'---.
    !             I       |       |      |                      \
    !              O      v       v      v                       \________.   downstream water
    !               N                                 baseflow             \--------
    !                     Groundwater Storage -- -- -- -- -- -- -- -->      \_______
    !
    !------------------------------------------------------------------------


    ! baseflow (m3) from the groundwater storage (m3) of the lake
    baseflow = baseflow_recession_coefficient * lake_gw_storage_previous   ![m3 in "hydrological" TS]

    ! Water balance for Groundwater storage
    ! percolation (mm to m3) added to groundwater storage (m3)
    lake_gw_storage = lake_gw_storage_previous + perc / 1000.0_dp * lake_surface_area - baseflow

    ! baseflow in m3/s
    baseflow = baseflow / tsHydroSecs


  END SUBROUTINE lake_gws_baseflow

END MODULE mo_lake_groundwaterstorage
