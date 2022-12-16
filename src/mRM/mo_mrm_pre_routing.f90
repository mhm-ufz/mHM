!> \file mo_mrm_pre_routing.f90
!> \brief \copybrief mo_mrm_pre_routing
!> \details \copydetails mo_mrm_pre_routing

!> \brief Performs pre-processing for routing for mHM at level L11.
!> \details This module performs runoff accumulation from L1 to L11 and inflow summation.
!> \changelog
!! - Stephan Thober Aug 2015
!!   - adapted to mRM
!! - Sebastian Mueller Jun 2020
!!   - separate module for pre-processing
!> \authors Luis Samaniego
!> \date Dec 2012
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
MODULE mo_mrm_pre_routing

  ! This module performs  pre-processing for routing for mHM at level L11.

  ! Written Sebastian Mueller Jun 2020

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: L11_runoff_acc
  PUBLIC :: add_inflow
  PUBLIC :: L11_meteo_acc
  PUBLIC :: calc_L1_runoff_E

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_runoff_acc

  !    PURPOSE
  !>       \brief total runoff accumulation at L11.

  !>       \details Upscales runoff in space from L1 to L11 if routing resolution
  !>       is higher than hydrology resolution (map_flag equals .true.) or
  !>       downscales runoff from L1 to L11 if routing resolution is lower
  !>       than hydrology resolution.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: qall"         total runoff L1 [mm TS-1]
  !>       \param[in] "real(dp), dimension(:) :: efecarea"     effective area in [km2] at Level 1
  !>       \param[in] "integer(i4), dimension(:) :: L1_L11_Id" L11 Ids mapped on L1
  !>       \param[in] "real(dp), dimension(:) :: L11_areacell" effective area in [km2] at Level 11
  !>       \param[in] "integer(i4), dimension(:) :: L11_L1_Id" L1 Ids mapped on L11
  !>       \param[in] "integer(i4) :: TS"                      time step in [s]
  !>       \param[in] "logical :: map_flag"                    Flag indicating whether routing resolution is higher than
  !>       hydrologic one

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: qAcc" aggregated runoff at L11 [m3 s-1]

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Jan 2013

  ! Modifications:
  ! Matthias Zink  Mar 2014 - added inflow from upstream areas
  ! Matthias Zink  Dec 2014 - adopted inflow gauges to ignore headwater cells
  ! Stephan Thober Sep 2015 - included downscaling of runoff
  ! Stephan Thober Feb 2016 - refactored upscaling of discharge from L1 to L11
  ! Stephan Thober Feb 2016 - refactored downscaling of discharge from L1 to L11
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  SUBROUTINE L11_runoff_acc(qAll, efecArea, L1_L11_Id, L11_areaCell, L11_L1_Id, TS, map_flag, qAcc)

    use mo_constants, only : HourSecs
    use mo_common_constants, only : nodata_dp

    implicit none

    ! total runoff L1 [mm TS-1]
    real(dp), intent(in), dimension(:) :: qall
    ! effective area in [km2] at Level 1
    real(dp), intent(in), dimension(:) :: efecarea
    ! L11 Ids mapped on L1
    integer(i4), intent(in), dimension(:) :: L1_L11_Id
    ! effective area in [km2] at Level 11
    real(dp), intent(in), dimension(:) :: L11_areacell
    ! L1 Ids mapped on L11
    integer(i4), intent(in), dimension(:) :: L11_L1_Id
    ! time step in [h]
    integer(i4), intent(in) :: TS
    ! Flag indicating whether routing resolution is higher than hydrologic one
    logical, intent(in) :: map_flag
    ! aggregated runoff at L11 [m3 s-1]
    real(dp), intent(out), dimension(:) :: qAcc

    integer(i4) :: k
    ! [s] time step
    real(dp) :: TST


    ! ------------------------------------------------------------------
    ! ACCUMULATION OF DISCHARGE TO A ROUTING CELL
    ! ------------------------------------------------------------------
    ! Hydrologic timestep in seconds
    TST = HourSecs * TS

    if (map_flag) then
      ! Estimate specific runoff at  L11
      ! NOTE:
      ! 1) Total discharge depth aggregated at L11 level [mm/TST]
      ! 2) Transform  depth [mm/TST] to discharge [m3/s]
      ! Total runoff should be divided by total_area to get
      ! specific discharge at L11. Then, to transform specific
      ! discharge from [mm/TST] to [m3/s], it should be multiplied by
      ! total_area [km2]*10^3 and divided by TST.
      ! Therefore, in this operation total_area cancels out.
      qAcc = 0._dp
      ! loop over high-resolution cells (L1) and add discharge to
      ! corresponding low-resolution cells (L11)
      do k = 1, size(qAll, 1)
        qAcc(L1_L11_Id(k)) = qAcc(L1_L11_Id(k)) + qAll(k) * efecArea(k)
      end do
      ! factor 1000 for conversion of mm*km2 in m3
      qAcc = qAcc * 1000.0_dp / TST
      !
    else
      ! initialize qout
      qAcc = nodata_dp
      do k = 1, size(qAcc, 1)
        ! map flux from coarse L1 resolution to fine L11 resolution
        qAcc(k) = qAll(L11_L1_Id(k))
      end do
      ! adjust flux by area cell
      ! factor 1000 for conversion of mm*km2 in m3
      qAcc(:) = qAcc(:) * L11_areaCell(:) * 1000.0_dp / TST
    end if

  END SUBROUTINE L11_runoff_acc

  ! ------------------------------------------------------------------

  !    NAME
  !        add_inflow

  !    PURPOSE
  !>       \brief
  !>          Adds inflow discharge to the runoff produced at the
  !>          cell where the inflow is occurring.
  !>       \details
  !>          If a inflow gauge is given, then this routine is adding the
  !>          values to the runoff produced at the grid cell where the
  !>          inflow is happening. The values are not directly added to the
  !>          river network. If this cell is not a headwater then the streamflow
  !>          produced upstream will be neglected.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: nInflowGauges"                 [-] number of inflow points
  !>       \param[in] "integer(i4), dimension(:) :: InflowIndexList" [-] index of inflow points
  !>       \param[in] "logical, dimension(:) :: InflowHeadwater"     [-] if to consider headwater cells of inflow gauge
  !>       \param[in] "integer(i4), dimension(:) :: InflowNodeList"  [-]        L11 ID of inflow points
  !>       \param[in] "real(dp), dimension(:) :: QInflow"            [m3 s-1]   inflowing water

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp), dimension(:) :: qOut" [m3 s-1] Series of attenuated runoff

  !    HISTORY
  !>       \authors Stephan Thober & Matthias Zink

  !>       \date Jul 2016

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine add_inflow(nInflowGauges, InflowIndexList, InflowHeadwater, InflowNodeList, QInflow, qOut)

    use mo_kind, only : dp, i4

    implicit none

    ! [-] number of inflow points
    integer(i4), intent(in) :: nInflowGauges
    ! [-] index of inflow points
    integer(i4), intent(in), dimension(:) :: InflowIndexList
    ! [-] if to consider headwater cells of inflow gauge
    logical, intent(in), dimension(:) :: InflowHeadwater
    ! [-]        L11 ID of inflow points
    integer(i4), intent(in), dimension(:) :: InflowNodeList
    ! [m3 s-1]   inflowing water
    real(dp), intent(in), dimension(:) :: QInflow
    ! [m3 s-1] Series of attenuated runoff
    real(dp), intent(inout), dimension(:) :: qOut

    integer(i4) :: ii


    ! discharge for inflow gauges (e.g. for missing upstream catchments) is added here
    ! should be put after UH attenuation because it is measured runoff at this cell
    if (nInflowGauges .gt. 0) then
      do ii = 1, nInflowGauges
        if (InflowHeadwater(ii)) then
          ! add inflowing water to water produced by upstream/headwater cells
          qOut(InflowNodeList(ii)) = qOut(InflowNodeList(ii)) + QInflow(InflowIndexList(ii))
        else
          ! put only timeseries and cut upstream/headwater cells produced water for routing
          qOut(InflowNodeList(ii)) = QInflow(InflowIndexList(ii))
        end if
      end do
    end if
  end subroutine add_inflow

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_E_acc

  !    PURPOSE
  !>       \brief temperature energy accumulation at L11.

  !>       \details Upscales energy in space from L1 to L11 if routing resolution
  !>       is higher than hydrology resolution (map_flag equals .true.) or
  !>       downscales runoff from L1 to L11 if routing resolution is lower
  !>       than hydrology resolution.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: qall"         total runoff L1 [mm K TS-1]
  !>       \param[in] "real(dp), dimension(:) :: efecarea"     effective area in [km2] at Level 1
  !>       \param[in] "integer(i4), dimension(:) :: L1_L11_Id" L11 Ids mapped on L1
  !>       \param[in] "real(dp), dimension(:) :: L11_areacell" effective area in [km2] at Level 11
  !>       \param[in] "integer(i4), dimension(:) :: L11_L1_Id" L1 Ids mapped on L11
  !>       \param[in] "integer(i4) :: TS"                      time step in [s]
  !>       \param[in] "logical :: map_flag"                    Flag indicating whether routing resolution is higher than
  !>       hydrologic one

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: qAcc" aggregated runoff at L11 [m3 K s-1]

  !    HISTORY
  !>       \authors Luis Samaniego

  !>       \date Jan 2013

  ! Modifications:
  ! Matthias Zink  Mar 2014 - added inflow from upstream areas
  ! Matthias Zink  Dec 2014 - adopted inflow gauges to ignore headwater cells
  ! Stephan Thober Sep 2015 - included downscaling of runoff
  ! Stephan Thober Feb 2016 - refactored upscaling of discharge from L1 to L11
  ! Stephan Thober Feb 2016 - refactored downscaling of discharge from L1 to L11
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  SUBROUTINE L11_E_acc(qAll, efecArea, L1_L11_Id, L11_areaCell, L11_L1_Id, TS, map_flag, qAcc)

    use mo_constants, only : HourSecs
    use mo_common_constants, only : nodata_dp

    implicit none

    ! total runoff L1 [mm TS-1]
    real(dp), intent(in), dimension(:) :: qall
    ! effective area in [km2] at Level 1
    real(dp), intent(in), dimension(:) :: efecarea
    ! L11 Ids mapped on L1
    integer(i4), intent(in), dimension(:) :: L1_L11_Id
    ! effective area in [km2] at Level 11
    real(dp), intent(in), dimension(:) :: L11_areacell
    ! L1 Ids mapped on L11
    integer(i4), intent(in), dimension(:) :: L11_L1_Id
    ! time step in [s]
    integer(i4), intent(in) :: TS
    ! Flag indicating whether routing resolution is higher than hydrologic one
    logical, intent(in) :: map_flag
    ! aggregated runoff at L11 [m3 s-1]
    real(dp), intent(out), dimension(:) :: qAcc

    integer(i4) :: k
    ! [s] time step
    real(dp) :: TST


    ! ------------------------------------------------------------------
    ! ACCUMULATION OF DISCHARGE TO A ROUTING CELL
    ! ------------------------------------------------------------------
    ! Hydrologic timestep in seconds
    TST = HourSecs * TS

    if (map_flag) then
      ! Estimate specific runoff at  L11
      ! NOTE:
      ! 1) Total discharge depth aggregated at L11 level [mm/TST]
      ! 2) Transform  depth [mm/TST] to discharge [m3/s]
      ! Total runoff should be divided by total_area to get
      ! specific discharge at L11. Then, to transform specific
      ! discharge from [mm/TST] to [m3/s], it should be multiplied by
      ! total_area [km2]*10^3 and divided by TST.
      ! Therefore, in this operation total_area cancels out.
      qAcc = 0._dp
      ! loop over high-resolution cells (L1) and add discharge to
      ! corresponding low-resolution cells (L11)
      do k = 1, size(qAll, 1)
        qAcc(L1_L11_Id(k)) = qAcc(L1_L11_Id(k)) + qAll(k) * efecArea(k)
      end do
      qAcc = qAcc * 1000.0_dp / TST
      !
    else
      ! initialize qout
      qAcc = nodata_dp
      do k = 1, size(qAcc, 1)
        ! map flux from coarse L1 resolution to fine L11 resolution
        qAcc(k) = qAll(L11_L1_Id(k))
      end do
      ! adjust flux by area cell
      qAcc(:) = qAcc(:) * L11_areaCell(:) * 1000.0_dp / TST
    end if

  END SUBROUTINE L11_E_acc

  ! ------------------------------------------------------------------

  !    NAME
  !        calc_L1_runoff_E

  !    PURPOSE
  !>       \brief calculate lateral temperature energy from runoff components.

  !>       \details calculate lateral temperature energy from runoff components.

  !    INTENT(IN)
  !>       \param[in] "REAL(dp) :: fSealed_area_fraction" sealed area fraction [1]
  !>       \param[in] "REAL(dp) :: fast_interflow"        \f$ q_0 \f$ Fast runoff component [mm TS-1]
  !>       \param[in] "REAL(dp) :: slow_interflow"        \f$ q_1 \f$ Slow runoff component [mm TS-1]
  !>       \param[in] "REAL(dp) :: baseflow"              \f$ q_2 \f$ Baseflow [mm TS-1]
  !>       \param[in] "REAL(dp) :: direct_runoff"         \f$ q_D \f$ Direct runoff from impervious areas  [mm TS-1]
  !>       \param[in] "REAL(dp) :: temp_air"              air temperature [K]
  !>       \param[in] "REAL(dp) :: mean_temp_air"         annual mean air temperature  [K]

  !    INTENT(OUT)
  !>       \param[out] "REAL(dp) :: lateral_E" \f$ E_T \f$ Generated runoff [K mm TS-1]

  !    HISTORY
  !>       \authors Sebastian Mueller

  !>       \date Jun 2020

  SUBROUTINE calc_L1_runoff_E( &
    fSealed_area_fraction, &
    fast_interflow, &
    slow_interflow, &
    baseflow, &
    direct_runoff, &
    temp_air, &
    mean_temp_air, &
    lateral_E &
  )

    use mo_constants, only : T0_dp  ! 273.15 - Celcius <-> Kelvin [K]

    implicit none

    ! sealed area fraction [1]
    REAL(dp), dimension(:), INTENT(IN) :: fSealed_area_fraction
    ! \f$ q_0 \f$ Fast runoff component [mm TS-1]
    REAL(dp), dimension(:), INTENT(IN) :: fast_interflow
    ! \f$ q_1 \f$ Slow runoff component [mm TS-1]
    REAL(dp), dimension(:), INTENT(IN) :: slow_interflow
    ! \f$ q_2 \f$ Baseflow [mm TS-1]
    REAL(dp), dimension(:), INTENT(IN) :: baseflow
    ! \f$ q_D \f$ Direct runoff from impervious areas  [mm TS-1]
    REAL(dp), dimension(:), INTENT(IN) :: direct_runoff
    ! air temperature [degC]
    real(dp), dimension(:), intent(in) :: temp_air
    ! annual mean air temperature [degC]
    real(dp), dimension(:), intent(in) :: mean_temp_air
    ! \f$ E_T \f$ Generated lateral Energy [K mm TS-1]
    REAL(dp), dimension(:), INTENT(inout) :: lateral_E

    ! convert temperatures form [deg C] to [K]
    ! accumulate in [K mm TS-1] -> convert later to [K m3 s-1] on L11
    ! following Wanders et.al. 2019
    lateral_E = lateral_E + ( &
      (baseflow * max(T0_dp + 5.0_dp, mean_temp_air + T0_dp) &
        + (slow_interflow + fast_interflow) * max(T0_dp, temp_air + T0_dp)) &
      * (1.0_dp - fSealed_area_fraction) &
      + direct_runoff * max(T0_dp, temp_air + T0_dp - 1.5_dp) * fSealed_area_fraction &
    )

  END SUBROUTINE calc_L1_runoff_E

  ! ------------------------------------------------------------------

  !    NAME
  !        L11_meteo_acc

  !    PURPOSE
  !>       \brief meteo forcing accumulation at L11 for temperature routing.

  !>       \details Upscales meteo forcing in space from L1 to L11 if routing resolution
  !>       is higher than hydrology resolution (map_flag equals .true.) or
  !>       downscales meteo forcing from L1 to L11 if routing resolution is lower
  !>       than hydrology resolution.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: meteo_all"    meteo forcing
  !>       \param[in] "real(dp), dimension(:) :: efecarea"     effective area in [km2] at Level 1
  !>       \param[in] "integer(i4), dimension(:) :: L1_L11_Id" L11 Ids mapped on L1
  !>       \param[in] "real(dp), dimension(:) :: L11_areacell" effective area in [km2] at Level 11
  !>       \param[in] "integer(i4), dimension(:) :: L11_L1_Id" L1 Ids mapped on L11
  !>       \param[in] "logical :: map_flag"                    Flag indicating whether routing resolution is higher than
  !>       hydrologic one

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: meteo_acc" aggregated meteo forcing

  !    HISTORY
  !>       \authors Sebastian Mueller

  !>       \date Jul 2020

  SUBROUTINE L11_meteo_acc(meteo_all, efecArea, L1_L11_Id, L11_areaCell, L11_L1_Id, map_flag, meteo_acc)

    use mo_common_constants, only : nodata_dp

    implicit none

    ! meteo forcing (state variable)
    real(dp), intent(in), dimension(:) :: meteo_all
    ! effective area in [km2] at Level 1
    real(dp), intent(in), dimension(:) :: efecarea
    ! L11 Ids mapped on L1
    integer(i4), intent(in), dimension(:) :: L1_L11_Id
    ! effective area in [km2] at Level 11
    real(dp), intent(in), dimension(:) :: L11_areacell
    ! L1 Ids mapped on L11
    integer(i4), intent(in), dimension(:) :: L11_L1_Id
    ! Flag indicating whether routing resolution is higher than hydrologic one
    logical, intent(in) :: map_flag
    ! aggregated meteo forcing
    real(dp), intent(out), dimension(:) :: meteo_acc

    integer(i4) :: k

    if (map_flag) then
      meteo_acc = 0._dp
      ! loop over high-resolution cells (L1) and add meteo forcing (weighted by area) to
      ! corresponding low-resolution cells (L11)
      do k = 1, size(meteo_all, 1)
        meteo_acc(L1_L11_Id(k)) = meteo_acc(L1_L11_Id(k)) + meteo_all(k) * efecArea(k)
      end do
      ! divide by L11 cell area to get weighted mean
      meteo_acc = meteo_acc / L11_areaCell(:)
    else
      ! initialize
      meteo_acc = nodata_dp
      do k = 1, size(meteo_acc, 1)
        ! map from coarse L1 resolution to fine L11 resolution
        meteo_acc(k) = meteo_all(L11_L1_Id(k))
      end do
    end if

  END SUBROUTINE L11_meteo_acc

END MODULE mo_mrm_pre_routing
