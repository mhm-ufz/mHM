!> \file mo_mlm_mpr.f90

!> \brief lake multiscale parameter regionalization module

!> \details  Detached lake multiscale parameter regionalization (MPR) for the
!>            saturated hydraulic conductivity of the lakebed. The lakebed Ks is
!>            regionalized per lake from the soil texture (sand and clay) over the
!>            lake footprint, using the Cosby et al. (1984) pedotransfer function,
!>            and then aggregated to a single value per lake by an arithmetic (cell-count)
!>            mean over the lake footprint cells. This mirrors the
!>            detached mRM regionalization (routing parameters + level-0 river
!>            slope -> celerity): it owns its own input (soil texture on the lake
!>            footprint) and target (one value per lake) and feeds a pedotransfer
!>            function, independent of the land-surface soil-moisture MPR.
!>            The resulting per-lake Ks replaces the former global lakebed Ks
!>            stand-in in the Darcy lakebed leakance Ks/b of the lake seepage.

!> \authors Juliane Mai, Stephan Thober, Rohini Kumar, Giovanni Dalmasso, Pallav K Shrestha, Sebastian Mueller

!> \date December 2012

! Modifications:
! Pallav K Shrestha, October 2026 - detached lake-MPR module (mo_mlm_mpr): the per-lake lakebed
!     conductivity regionalizer lake_derive_ks reuses the per-lake footprint aggregation of
!     mo_lake_topology:lake_derive_max_surface_area_dem (Pallav Shrestha, Sebastian Mueller) and the
!     Cosby (1984) saturated-conductivity pedotransfer function of mo_mpr_soilmoist:hydro_cond
!     (Stephan Thober, Rohini Kumar); the soil look-up table reader lake_read_soil_texture_lut reuses
!     the parse logic of mo_soil_database:read_soil_LUT (Juliane Mai); the per-lake upscaling
!     upscale_arithmetic_mean_lake mirrors the level-0 arithmetic (cell-count) mean of
!     mo_upscaling_operators:upscale_arithmetic_mean (Giovanni Dalmasso, Rohini Kumar), evaluated
!     over the packed lake footprint instead of a level-1 cell window. Assembled for the mLM v6
!     migration so the lake conductivity is regionalized independently of the land-surface soil-moisture MPR.

#include "logging.h"

MODULE mo_mlm_mpr

  USE mo_logging
  USE mo_kind,   only : i4, i8, dp
  USE mo_grid,   only : grid_t
  IMPLICIT NONE

  PUBLIC :: lake_derive_ks            ! regionalizes the per-lake lakebed Ks from soil texture over the lake footprint
  PUBLIC :: lake_read_soil_texture_lut ! reads top-horizon sand/clay per soil class from the soil LUT
  PUBLIC :: upscale_arithmetic_mean_lake ! per-lake arithmetic (cell-count) mean of a level-0 lake-footprint field

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         lake_derive_ks

  !     PURPOSE
  !>        \brief Derives the per-lake lakebed saturated hydraulic conductivity.

  !>        \details Computes the lakebed saturated hydraulic conductivity Ks for
  !>                 every lake cell over the lake footprint with the Cosby et al.
  !>                 (1984) pedotransfer function from the top-horizon sand and clay
  !>                 content, then aggregates it to one value per lake by an
  !>                 arithmetic (cell-count) mean over the footprint cells (the same operator the
  !>                 land-surface MPR uses, upscale_arithmetic_mean). The
  !>                 footprint cells are labelled by their stable lake ID in lake_map
  !>                 (aligned with the packed level-0 lake grid) and matched to the
  !>                 lake points through lake_ids, exactly as in the DEM-footprint
  !>                 area aggregation (lake_derive_max_surface_area_dem). The per-cell
  !>                 sand and clay fractions are supplied by the caller (resolved from the
  !>                 soil look-up table over the lake footprint by the Input container), so this
  !>                 routine owns no global connections and is self-contained.

  !     CALLING SEQUENCE
  !         call lake_derive_ks(grid, lake_map, lake_ids, sand, clay, ptf_ks, ks_lake)

  !     INTENT(IN)
  !>       \param[in] "type(grid_t)             :: grid"      packed level-0 lake footprint grid (carries cell_area)
  !>       \param[in] "integer(i8), dimension(:) :: lake_map"  [1]  stable lake ID on each packed footprint cell
  !>       \param[in] "integer(i8), dimension(:) :: lake_ids"  [1]  stable lake ID in lake-point order
  !>       \param[in] "real(dp),    dimension(:) :: sand"      [%]  top-horizon sand content on each footprint cell
  !>       \param[in] "real(dp),    dimension(:) :: clay"      [%]  top-horizon clay content on each footprint cell
  !>       \param[in] "real(dp),    dimension(4) :: ptf_ks"         Cosby Ks pedotransfer parameters:
  !>                                                                 (1) constant, (2) sand coefficient,
  !>                                                                 (3) clay coefficient, (4) unit conversion (inch/h -> cm/d)

  !     INTENT(OUT)
  !>       \param[out] "real(dp), allocatable, dimension(:) :: ks_lake"  [cm d-1] arithmetic-mean lakebed Ks per lake

  !     RESTRICTIONS
  !>       \note The Cosby pedotransfer function and its minimum floor (1.1 cm/d) are
  !>             kept identical to the land-surface implementation (mo_mpr_soilmoist:hydro_cond)
  !>             so the lakebed and land-surface Ks share the same definition.

  !     LITERATURE
  !>       Cosby, B. J., Hornberger, G. M., Clapp, R. B., & Ginn, T. R. (1984).
  !>       A statistical exploration of the relationships of soil moisture
  !>       characteristics to the physical properties of soils. WRR, 20(6), 682-690.

  !     HISTORY
  !>       \authors Stephan Thober, Rohini Kumar (Cosby Ks pedotransfer function, from mo_mpr_soilmoist)
  !>       \date    December 2012
  !        Modifications:
  !        Pallav K Shrestha, October 2026 - applied the Cosby pedotransfer function per lake-footprint
  !            cell and upscaled it per lake with upscale_arithmetic_mean_lake, for the mLM v6 migration.

  subroutine lake_derive_ks(grid, lake_map, lake_ids, sand, clay, ptf_ks, ks_lake)

    implicit none

    type(grid_t),              intent(in)  :: grid        ! packed level-0 lake footprint grid
    integer(i8), dimension(:), intent(in)  :: lake_map    ! stable lake ID on each footprint cell
    integer(i8), dimension(:), intent(in)  :: lake_ids    ! stable lake ID in lake-point order
    real(dp),    dimension(:), intent(in)  :: sand        ! [%] sand content on each footprint cell (top horizon)
    real(dp),    dimension(:), intent(in)  :: clay        ! [%] clay content on each footprint cell (top horizon)
    real(dp),    dimension(4), intent(in)  :: ptf_ks      ! Cosby Ks pedotransfer parameters
    real(dp), allocatable, dimension(:), intent(out) :: ks_lake  ! [cm d-1] per-lake mean lakebed Ks

    ! local variables
    real(dp), allocatable :: ks_cell(:)    ! [cm d-1] Cosby Ks on each footprint cell
    real(dp) :: x                          ! temporary Cosby exponent argument
    real(dp), parameter :: Ks_c = 10.0_dp  ! Cosby base (matches mo_mpr_constants)
    real(dp), parameter :: Ks_floor = 1.1_dp ! [cm d-1] minimum Ks (matches mo_mpr_soilmoist:hydro_cond)
    integer(i8) :: i

    if (size(sand, kind=i8) /= grid%ncells .or. size(clay, kind=i8) /= grid%ncells) then
      log_fatal(*) "Lake MPR: soil texture size does not match the lake footprint grid."
      error stop 1
    end if

    ! Cosby et al. (1984) saturated hydraulic conductivity [cm d-1] per footprint cell, identical to the
    ! land-surface pedotransfer function (mo_mpr_soilmoist:hydro_cond):
    !   Ks = ptf(4) * exp( (ptf(1) + ptf(2)*sand - ptf(3)*clay) * ln(Ks_c) )
    allocate(ks_cell(grid%ncells))
    do i = 1_i8, grid%ncells
      x = ptf_ks(1) + ptf_ks(2) * sand(i) - ptf_ks(3) * clay(i)
      ks_cell(i) = ptf_ks(4) * exp(x * log(Ks_c))
      if (ks_cell(i) < Ks_floor) ks_cell(i) = Ks_floor
    end do

    ! Upscale the level-0 lakebed Ks to one value per lake by the arithmetic (cell-count) mean over
    ! the lake footprint, the same operator the land-surface MPR uses (upscale_arithmetic_mean).
    call upscale_arithmetic_mean_lake(lake_map, lake_ids, ks_cell, ks_lake)

    deallocate(ks_cell)

  end subroutine lake_derive_ks


  ! ------------------------------------------------------------------

  !     NAME
  !         lake_read_soil_texture_lut

  !     PURPOSE
  !>        \brief Reads the top-horizon sand and clay content per soil class from the soil LUT.

  !>        \details Parses the classical mHM soil class-definition look-up table and returns, for
  !>                 every soil class, the sand and clay content of its uppermost horizon (horizon 1).
  !>                 The file format is: a first line "nSoil_Types <n>", a header line, then one row
  !>                 per class-and-horizon "<class> <horizon> <upper_depth> <lower_depth> <clay%> <sand%>
  !>                 <bulk_density>". Only horizon 1 of each class is retained (the lakebed interfaces
  !>                 the surface soil). Classes absent from the table keep the missing-value sentinel.
  !>                 This is a self-contained reader so the detached lake MPR does not depend on the
  !>                 land-surface soil database.

  !     CALLING SEQUENCE
  !         call lake_read_soil_texture_lut(path, n_classes, sand_by_class, clay_by_class)

  !     INTENT(IN)
  !>       \param[in] "character(*) :: path"          soil LUT file path

  !     INTENT(OUT)
  !>       \param[out] "integer(i4) :: n_classes"     number of soil classes declared in the LUT header
  !>       \param[out] "real(dp), allocatable :: sand_by_class(:)" [%] horizon-1 sand content per class (missing = nodata)
  !>       \param[out] "real(dp), allocatable :: clay_by_class(:)" [%] horizon-1 clay content per class (missing = nodata)

  !     HISTORY
  !>       \authors Juliane Mai (soil look-up table parsing, from mo_soil_database:read_soil_LUT)
  !>       \date    December 2012
  !        Modifications:
  !        Pallav K Shrestha, October 2026 - reused the soil-LUT parse logic as a self-contained reader
  !            returning only the top-horizon sand/clay per class, for the detached lake MPR.

  subroutine lake_read_soil_texture_lut(path, n_classes, sand_by_class, clay_by_class)

    use mo_constants, only : nodata_dp

    implicit none

    character(*),              intent(in)  :: path           ! soil LUT path
    integer(i4),               intent(out) :: n_classes      ! number of soil classes
    real(dp), allocatable,     intent(out) :: sand_by_class(:) ! [%] horizon-1 sand per class
    real(dp), allocatable,     intent(out) :: clay_by_class(:) ! [%] horizon-1 clay per class

    integer :: unit               ! default-kind unit for newunit=
    integer(i4) :: ios
    integer(i4) :: id, horizon
    real(dp)    :: up, down, cly, snd, bd
    character(256) :: dummy
    logical :: exists

    inquire(file=path, exist=exists)
    if (.not.exists) then
      log_fatal(*) "Lake MPR: soil LUT not found: ", trim(path)
      error stop 1
    end if

    open(newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0_i4) then
      log_fatal(*) "Lake MPR: cannot open soil LUT: ", trim(path)
      error stop 1
    end if

    ! header line: "nSoil_Types <n>"
    read(unit, *, iostat=ios) dummy, n_classes
    if (ios /= 0_i4 .or. n_classes < 1_i4) then
      log_fatal(*) "Lake MPR: malformed soil LUT header (expected 'nSoil_Types <n>'): ", trim(path)
      error stop 1
    end if
    ! column header line
    read(unit, *, iostat=ios) dummy

    allocate(sand_by_class(n_classes), source=nodata_dp)
    allocate(clay_by_class(n_classes), source=nodata_dp)

    ! data rows: <class> <horizon> <up> <down> <clay> <sand> <bd>
    do
      read(unit, *, iostat=ios) id, horizon, up, down, cly, snd, bd
      if (ios /= 0_i4) exit
      if (id < 1_i4 .or. id > n_classes) then
        log_fatal(*) "Lake MPR: soil LUT class id out of range in ", trim(path)
        error stop 1
      end if
      ! keep the uppermost horizon only
      if (horizon == 1_i4) then
        sand_by_class(id) = snd
        clay_by_class(id) = cly
      end if
    end do

    close(unit)

  end subroutine lake_read_soil_texture_lut


  ! ------------------------------------------------------------------

  !     NAME
  !         upscale_arithmetic_mean_lake

  !     PURPOSE
  !>        \brief Per-lake arithmetic (cell-count) mean of a level-0 lake-footprint field.

  !>        \details Upscales a level-0 field defined on the packed lake footprint to one value per
  !>                 lake by the arithmetic mean, giving every footprint cell equal weight:
  !>                 ks_lake(k) = sum( field over the cells of lake k ) / (number of those cells).
  !>                 Footprint cells are labelled by their stable lake ID in lake_map and matched to
  !>                 the lake points through lake_ids. This is the lake-footprint counterpart of the
  !>                 land-surface level-0 -> level-1 operator mo_upscaling_operators:upscale_arithmetic_mean:
  !>                 the same cell-count arithmetic mean, evaluated over a lake footprint rather than a
  !>                 level-1 cell window (on a regular grid the two are identical).

  !     CALLING SEQUENCE
  !         call upscale_arithmetic_mean_lake(lake_map, lake_ids, l0_field, lake_mean)

  !     INTENT(IN)
  !>       \param[in] "integer(i8), dimension(:) :: lake_map"  [1] stable lake ID on each packed footprint cell
  !>       \param[in] "integer(i8), dimension(:) :: lake_ids"  [1] stable lake ID in lake-point order
  !>       \param[in] "real(dp),    dimension(:) :: l0_field"      level-0 field on the lake footprint cells

  !     INTENT(OUT)
  !>       \param[out] "real(dp), allocatable, dimension(:) :: lake_mean"  per-lake arithmetic mean

  !     HISTORY
  !>       \authors Giovanni Dalmasso, Rohini Kumar (arithmetic cell-count mean, from mo_upscaling_operators)
  !>       \date    December 2012
  !        Modifications:
  !        Pallav K Shrestha, September 2018 - evaluated over the packed lake footprint (lake_map) rather than
  !            a level-1 cell window, for the detached lake MPR.

  subroutine upscale_arithmetic_mean_lake(lake_map, lake_ids, l0_field, lake_mean)

    implicit none

    integer(i8), dimension(:), intent(in)  :: lake_map   ! stable lake ID on each footprint cell
    integer(i8), dimension(:), intent(in)  :: lake_ids   ! stable lake ID in lake-point order
    real(dp),    dimension(:), intent(in)  :: l0_field   ! level-0 field on the lake footprint cells
    real(dp), allocatable, dimension(:), intent(out) :: lake_mean  ! per-lake arithmetic mean

    ! local variables
    integer(i8), allocatable :: cell_count(:)  ! number of footprint cells per lake
    integer(i8) :: i, j

    if (size(l0_field, kind=i8) /= size(lake_map, kind=i8)) then
      log_fatal(*) "Lake MPR: field size does not match the lake map."
      error stop 1
    end if

    allocate(lake_mean(size(lake_ids)), source=0.0_dp)
    allocate(cell_count(size(lake_ids)), source=0_i8)

    do i = 1_i8, size(lake_map, kind=i8)
      j = findloc(lake_ids, lake_map(i), dim=1, kind=i8)
      if (j < 1_i8) then
        log_fatal(*) "Lake MPR: lake map contains an unknown stable lake ID."
        error stop 1
      end if
      lake_mean(j)  = lake_mean(j)  + l0_field(i)
      cell_count(j) = cell_count(j) + 1_i8
    end do

    if (any(cell_count < 1_i8)) then
      log_fatal(*) "Lake MPR: every lake must have at least one footprint cell."
      error stop 1
    end if

    lake_mean = lake_mean / real(cell_count, dp)

    deallocate(cell_count)

  end subroutine upscale_arithmetic_mean_lake

END MODULE mo_mlm_mpr
