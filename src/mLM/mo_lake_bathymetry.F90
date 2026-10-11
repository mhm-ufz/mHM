!> \file mo_lake_bathymetry.f90

!> \brief lake bathymetry generation module

!> \details  Generates the height-area-volume (hAV) curve of a lake from a built-in
!>           analytical bathymetry shape, using the spill reference level, surface area
!>           and volume. Ported from mo_mlm_prepare :: prepare_lake_bathymetry (mLM v5.13).

!> \authors Pallav Kumar Shrestha
!> \date Jan 2023

MODULE mo_lake_bathymetry

  USE mo_kind,   only : i4, i8, dp

  IMPLICIT NONE

    PUBLIC :: generate_lake_bathymetry
    PUBLIC :: lake_hav_t

  ! ------------------------------------------------------------------

  !> \brief Per-lake height-area-volume (hAV) curve, tabulated at a 1 cm depth interval.
  !> \details One instance per lake; the three columns share a common length (the number
  !>          of depth intervals from bottom to crest, which differs between lakes). All in SI.
  type :: lake_hav_t
    real(dp), allocatable :: level(:)         !< hAV level column [masl]
    real(dp), allocatable :: surface_area(:)  !< hAV surface area column [m2]
    real(dp), allocatable :: volume(:)        !< hAV volume column [m3]
  end type lake_hav_t

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        generate_lake_bathymetry

  !    PURPOSE
  !>       \brief generate the hAV (height-area-volume) curve of a lake

  !>       \details Builds the lake bathymetry (hAV curve) from a built-in analytical
  !>                shape using the spill reference level, surface area and volume.
  !>                The curve is tabulated at a 1 cm depth interval from the lake bottom
  !>                up to the crest. All quantities are in SI units (m, m2, m3).
  !>                Shape types:
  !>                  1 - linear h~V relation, with constant area as the maximum area
  !>                  2 - Liebe (2005) pyramid with a square base diagonally cut into halves

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: bathymetry_type"      [-]      built-in shape selector
  !>       \param[in] "real(dp)    :: lake_spill_reference_level"       [masl]   spill reference level
  !>       \param[in] "real(dp)    :: lake_spill_reference_surface_area"[m2]     surface area at spill reference level
  !>       \param[in] "real(dp)    :: lake_spill_reference_volume"      [m3]     volume at spill reference level

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>       \param[out] "real(dp), allocatable :: hav_lake_level(:)"        [masl] hAV level column
  !>       \param[out] "real(dp), allocatable :: hav_lake_surface_area(:)" [m2]   hAV surface area column
  !>       \param[out] "real(dp), allocatable :: hav_lake_volume(:)"       [m3]   hAV volume column

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
  !         Liebe, J. (2005). Estimation of water storage capacity and evaporation losses
  !         of small reservoirs in the Upper East Region of Ghana.

  !     HISTORY
  !>        \author Pallav Shrestha
  !>        \date Jan 2023

  !     MODIFICATIONS
  !        Pallav Shrestha, Oct 2026 - mLM v6 migration: extracted from
  !                                    mo_mlm_prepare::prepare_lake_bathymetry into its own
  !                                    module. Converted to SI units (dropped the km2/mcm
  !                                    1e6 factors), consistent nomenclature, and passed
  !                                    inputs/outputs as arguments instead of the lake
  !                                    global structure. Corrected the final partial row so
  !                                    the crest row equals the reported maxima exactly.

  SUBROUTINE generate_lake_bathymetry (bathymetry_type, lake_spill_reference_level, lake_spill_reference_surface_area, &
                                       lake_spill_reference_volume, hav_lake_level, hav_lake_surface_area, &
                                       hav_lake_volume)

    USE mo_constants,             ONLY : ONETHIRD_dp
    USE mo_string_utils,          ONLY : num2str
    USE mo_message,               ONLY : message

    IMPLICIT NONE

    integer(i4),         intent(in)  :: bathymetry_type         ! built-in shape selector
    real(dp),            intent(in)  :: lake_spill_reference_level          ! spill reference level [masl]
    real(dp),            intent(in)  :: lake_spill_reference_surface_area   ! surface area at spill reference level [m2]
    real(dp),            intent(in)  :: lake_spill_reference_volume         ! volume at spill reference level [m3]
    real(dp), allocatable, intent(out) :: hav_lake_level(:)         ! hAV level column [masl]
    real(dp), allocatable, intent(out) :: hav_lake_surface_area(:)  ! hAV surface area column [m2]
    real(dp), allocatable, intent(out) :: hav_lake_volume(:)        ! hAV volume column [m3]

    ! local variables

    real(dp)          :: depth_max, f
    integer(i8)       :: nIntervals, nRows, i
    logical           :: has_partial_row

    ! Check reported inputs
    if (lake_spill_reference_surface_area <= 0._dp .or. lake_spill_reference_volume <= 0._dp) then
      call message('***ERROR: surface area at spill reference level and volume must be positive.')
      call message('          Error occured in mo_lake_bathymetry')
      stop 1
    end if

    ! Calculate the max depth (m) from the reported volume and surface area
    select case (bathymetry_type)
    case(1) ! linear h~V relation, with constant A as Amax
      depth_max = lake_spill_reference_volume / lake_spill_reference_surface_area
    case(2) ! Liebe (2005) pyramid with a square base diagonally cut into halves
      depth_max = 3._dp * lake_spill_reference_volume / lake_spill_reference_surface_area
    case default
      call message('***ERROR: unsupported bathymetry_type ', trim(adjustl(num2str(bathymetry_type))), '.')
      call message('          Error occured in mo_lake_bathymetry')
      stop 1
    end select

    ! Loop at 1 cm interval
    nIntervals = FLOOR(depth_max * 100)

    ! final partial row needed when depth_max is not an exact multiple of 1 cm
    has_partial_row = ( depth_max /= real(nIntervals, dp) / 100._dp )
    nRows = nIntervals + 1
    if (has_partial_row) nRows = nRows + 1

    allocate(hav_lake_level(nRows), hav_lake_surface_area(nRows), hav_lake_volume(nRows))

    select case (bathymetry_type)
    case(1) ! linear h~V relation, with constant A as Amax

      ! First entry
      hav_lake_level(1)        = 0._dp
      hav_lake_surface_area(:) = lake_spill_reference_surface_area
      hav_lake_volume(1)       = 0._dp

      do i = 2, nIntervals + 1
        hav_lake_level(i)  = hav_lake_level(i - 1) + 0.01_dp
        hav_lake_volume(i) = lake_spill_reference_surface_area * hav_lake_level(i)
      end do

      ! add the last entry at the crest
      if (has_partial_row) then
        hav_lake_level(nRows)  = depth_max
        hav_lake_volume(nRows) = lake_spill_reference_surface_area * depth_max
      end if

    case(2) ! Liebe (2005) pyramid with a square base diagonally cut into halves

      ! f - the ratio of characteristic length to depth
      ! Liebe (2005) proposes f = 20, which is not suitable for large reservoirs
      ! as the depth can reach 1 km! So get the value of f rather based on
      ! lake_spill_reference_volume, lake_spill_reference_surface_area, depth_max (SI: no km2 scaling)
      f = SQRT((lake_spill_reference_surface_area / depth_max ** 2._dp) ** 2._dp + 1._dp)

      ! First entry
      hav_lake_level(1)        = 0._dp
      hav_lake_surface_area(1) = 0._dp
      hav_lake_volume(1)       = 0._dp

      do i = 2, nIntervals + 1
        hav_lake_level(i)        = hav_lake_level(i - 1) + 0.01_dp
        hav_lake_surface_area(i) = hav_lake_level(i) ** 2._dp * SQRT( f**2._dp - 1 )
        hav_lake_volume(i)       = ONETHIRD_dp * hav_lake_level(i) ** 3._dp * SQRT( f**2._dp - 1 )
      end do

      ! add the last entry at the crest
      if (has_partial_row) then
        hav_lake_level(nRows)        = depth_max
        hav_lake_surface_area(nRows) = depth_max ** 2._dp * SQRT( f**2._dp - 1 )
        hav_lake_volume(nRows)       = ONETHIRD_dp * depth_max ** 3._dp * SQRT( f**2._dp - 1 )
      end if

    end select

    ! convert depth to elevation, placing the crest at the spill reference level
    hav_lake_level(:) = hav_lake_level(:) + lake_spill_reference_level - depth_max

  END SUBROUTINE generate_lake_bathymetry

  ! ------------------------------------------------------------------

END MODULE mo_lake_bathymetry
