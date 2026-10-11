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
  USE mo_string_utils, only : n2s => num2str
  IMPLICIT NONE

  PUBLIC :: lake_derive_ks            ! regionalizes the per-lake lakebed Ks from soil texture over the lake footprint
  PUBLIC :: lake_derive_kbf           ! regionalizes the per-lake baseflow recession coefficient from geology over the lake footprint
  PUBLIC :: lake_read_soil_texture_lut ! reads top-horizon sand/clay per soil class from the soil LUT
  PUBLIC :: lake_read_geology_lut     ! reads the geological unit IDs and their global baseflow-parameter index from the geology LUT
  PUBLIC :: lake_check_geo_units      ! checks that every geological unit on the lake footprint is present in the geology LUT
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
  !         lake_derive_kbf

  !     PURPOSE
  !>        \brief Derives the per-lake baseflow recession coefficient.

  !>        \details Regionalizes the lake groundwater-store baseflow recession coefficient k_bf for
  !>                 every lake from the geological units over the lake footprint, then aggregates it to
  !>                 one value per lake by an arithmetic (cell-count) mean over the footprint cells (the
  !>                 same operator the lakebed Ks uses, upscale_arithmetic_mean_lake / the land-surface
  !>                 upscale_arithmetic_mean). The regionalization is the geological baseflow-parameter
  !>                 lookup: each footprint cell carries a geological unit ID (geo_unit); the unit is
  !>                 located in the geology look-up table (geo_unit_list) and its global baseflow-parameter
  !>                 slot (geo_param_index) selects the per-geological-unit baseflow recession value from
  !>                 baseflow_recession. This is the v5 mpr_lake baseflow path (baseflow_param over the
  !>                 lake footprint + upscale_arithmetic_mean_lake) rebuilt in the detached lake-MPR,
  !>                 and matches the land-surface mpr_bridge_baseflow_param cell lookup exactly. The
  !>                 baseflow_recession values are supplied by the caller (the lake module's own lakes_1
  !>                 parameters, or the land baseflow process borrowed when it is active), so this routine
  !>                 owns no global connections and is self-contained.

  !     CALLING SEQUENCE
  !         call lake_derive_kbf(lake_map, lake_ids, geo_unit, geo_unit_list, geo_param_index, baseflow_recession, kbf_lake)

  !     INTENT(IN)
  !>       \param[in] "integer(i8), dimension(:) :: lake_map"           [1] stable lake ID on each packed footprint cell
  !>       \param[in] "integer(i8), dimension(:) :: lake_ids"           [1] stable lake ID in lake-point order
  !>       \param[in] "integer(i4), dimension(:) :: geo_unit"           [1] geological unit ID on each footprint cell
  !>       \param[in] "integer(i4), dimension(:) :: geo_unit_list"      [1] geological unit IDs in the geology LUT
  !>       \param[in] "integer(i4), dimension(:) :: geo_param_index"    [1] global baseflow-parameter slot per geology LUT class
  !>       \param[in] "real(dp),    dimension(:) :: baseflow_recession" [d] per-geological-unit baseflow recession values

  !     INTENT(OUT)
  !>       \param[out] "real(dp), allocatable, dimension(:) :: kbf_lake"  [d] arithmetic-mean baseflow recession coefficient per lake

  !     HISTORY
  !>       \authors Pallav Kumar Shrestha (lake baseflow-parameter regionalization, from v5 mo_multi_param_reg:mpr_lake / baseflow_param)
  !>       \date    October 2018
  !        Modifications:
  !        Pallav K Shrestha, October 2026 - rebuilt the v5 lake baseflow-parameter path in the detached
  !            lake-MPR: the per-footprint-cell geological-unit lookup (geo_unit -> geo_unit_list ->
  !            geo_param_index -> baseflow_recession) mirrors the land-surface mo_mpr_legacy_bridge:
  !            mpr_bridge_baseflow_param, upscaled per lake with upscale_arithmetic_mean_lake, for the
  !            mLM v6 migration.

  subroutine lake_derive_kbf(lake_map, lake_ids, geo_unit, geo_unit_list, geo_param_index, &
      baseflow_recession, kbf_lake)

    implicit none

    integer(i8), dimension(:), intent(in)  :: lake_map           ! stable lake ID on each footprint cell
    integer(i8), dimension(:), intent(in)  :: lake_ids           ! stable lake ID in lake-point order
    integer(i4), dimension(:), intent(in)  :: geo_unit           ! geological unit ID on each footprint cell
    integer(i4), dimension(:), intent(in)  :: geo_unit_list      ! geological unit IDs in the geology LUT
    integer(i4), dimension(:), intent(in)  :: geo_param_index    ! global baseflow-parameter slot per LUT class
    real(dp),    dimension(:), intent(in)  :: baseflow_recession ! [d] per-geological-unit baseflow recession
    real(dp), allocatable, dimension(:), intent(out) :: kbf_lake ! [d] per-lake mean baseflow recession coefficient

    ! local variables
    real(dp), allocatable :: kbf_cell(:)   ! [d] baseflow recession on each footprint cell
    integer(i4) :: class_pos               ! position of a cell's geo unit in geo_unit_list
    integer(i8) :: i
    integer(i4) :: j

    if (size(geo_param_index) /= size(geo_unit_list)) then
      log_fatal(*) "Lake MPR: geo_param_index and geology LUT unit list sizes do not match."
      error stop 1
    end if
    if (any(geo_param_index < 1_i4) .or. any(geo_param_index > size(baseflow_recession))) then
      log_fatal(*) "Lake MPR: geology LUT baseflow-parameter index is outside the baseflow parameters."
      error stop 1
    end if
    if (size(geo_unit, kind=i8) /= size(lake_map, kind=i8)) then
      log_fatal(*) "Lake MPR: geological unit field size does not match the lake footprint grid."
      error stop 1
    end if

    ! Per-footprint-cell baseflow recession coefficient, the geological baseflow-parameter lookup
    ! identical to the land-surface implementation (mo_mpr_legacy_bridge:mpr_bridge_baseflow_param):
    ! locate the cell's geological unit in the LUT unit list, then pick its global baseflow-parameter slot.
    allocate(kbf_cell(size(lake_map, kind=i8)))
    do i = 1_i8, size(lake_map, kind=i8)
      class_pos = 0_i4
      do j = 1_i4, size(geo_unit_list)
        if (geo_unit_list(j) == geo_unit(i)) then
          class_pos = j
          exit
        end if
      end do
      if (class_pos < 1_i4) then
        log_fatal(*) "Lake MPR: lake footprint geological unit ", n2s(geo_unit(i)), &
          " is not present in the geology LUT."
        error stop 1
      end if
      kbf_cell(i) = baseflow_recession(geo_param_index(class_pos))
    end do

    ! Upscale the level-0 baseflow recession to one value per lake by the arithmetic (cell-count) mean
    ! over the lake footprint, the same operator the lakebed Ks regionalization uses.
    call upscale_arithmetic_mean_lake(lake_map, lake_ids, kbf_cell, kbf_lake)

    deallocate(kbf_cell)

  end subroutine lake_derive_kbf


  ! ------------------------------------------------------------------

  !     NAME
  !         lake_read_geology_lut

  !     PURPOSE
  !>        \brief Reads the geological unit IDs and their global baseflow-parameter index from the geology LUT.

  !>        \details Parses the classical mHM geology class-definition look-up table and returns, for every
  !>                 geological class, its unit ID (the value stored in the geology raster) and its global
  !>                 baseflow-parameter index (the slot into the baseflow recession parameter set). The file
  !>                 format is: a first line "nGeo_Formations <n>", a header line, then one row per class
  !>                 "<GeoParam_index> <ClassUnit> <Karstic> <Description>". This is a self-contained reader
  !>                 (the parse logic of mo_read_lut:read_geoformation_lut) so the detached lake MPR does not
  !>                 depend on the land-surface geology class definition read by the MPR container; the lake
  !>                 module reads its own geology LUT and stays independent of the land-surface MPR.

  !     CALLING SEQUENCE
  !         call lake_read_geology_lut(path, n_classes, geo_unit_list, geo_param_index)

  !     INTENT(IN)
  !>       \param[in] "character(*) :: path"          geology LUT file path

  !     INTENT(OUT)
  !>       \param[out] "integer(i4) :: n_classes"                      number of geological classes in the LUT
  !>       \param[out] "integer(i4), allocatable :: geo_unit_list(:)"  geological unit IDs per class
  !>       \param[out] "integer(i4), allocatable :: geo_param_index(:)" global baseflow-parameter slot per class

  !     HISTORY
  !>       \authors Matthias Zink, Rohini Kumar (geology look-up table parsing, from mo_read_lut:read_geoformation_lut)
  !>       \date    January 2013
  !        Modifications:
  !        Pallav K Shrestha, October 2026 - reused the geology-LUT parse logic as a self-contained reader
  !            returning the unit IDs and the global baseflow-parameter index per class, for the detached lake MPR.

  subroutine lake_read_geology_lut(path, n_classes, geo_unit_list, geo_param_index)

    implicit none

    character(*),              intent(in)  :: path              ! geology LUT path
    integer(i4),               intent(out) :: n_classes         ! number of geological classes
    integer(i4), allocatable,  intent(out) :: geo_unit_list(:)  ! geological unit IDs per class
    integer(i4), allocatable,  intent(out) :: geo_param_index(:) ! global baseflow-parameter slot per class

    integer :: unit               ! default-kind unit for newunit=
    integer(i4) :: ios
    integer(i4) :: i, param_index, unit_id, karstic
    character(256) :: dummy
    logical :: exists

    inquire(file=path, exist=exists)
    if (.not.exists) then
      log_fatal(*) "Lake MPR: geology LUT not found: ", trim(path)
      error stop 1
    end if

    open(newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0_i4) then
      log_fatal(*) "Lake MPR: cannot open geology LUT: ", trim(path)
      error stop 1
    end if

    ! header line: "nGeo_Formations <n>"
    read(unit, *, iostat=ios) dummy, n_classes
    if (ios /= 0_i4 .or. n_classes < 1_i4) then
      log_fatal(*) "Lake MPR: malformed geology LUT header (expected 'nGeo_Formations <n>'): ", trim(path)
      error stop 1
    end if
    ! column header line
    read(unit, *, iostat=ios) dummy

    allocate(geo_unit_list(n_classes))
    allocate(geo_param_index(n_classes))

    ! data rows: <GeoParam_index> <ClassUnit> <Karstic> <Description>
    do i = 1_i4, n_classes
      read(unit, *, iostat=ios) param_index, unit_id, karstic, dummy
      if (ios /= 0_i4) then
        log_fatal(*) "Lake MPR: nGeo_Formations in the geology LUT exceeds the available rows: ", trim(path)
        error stop 1
      end if
      geo_param_index(i) = param_index
      geo_unit_list(i)   = unit_id
    end do

    close(unit)

  end subroutine lake_read_geology_lut


  ! ------------------------------------------------------------------

  !     NAME
  !         lake_check_geo_units

  !     PURPOSE
  !>        \brief Checks that every geological unit on the lake footprint is present in the geology LUT.

  !>        \details Mirrors the land-surface check (mo_mpr_container:mpr_check_geo_units_against_lut),
  !>                 but over the LAKE footprint rather than the land grid. The land-surface check only
  !>                 validates the geological units present on the level-0 land grid, which holes out the
  !>                 lake cells; a geological unit that occurs only under a lake is therefore never
  !>                 validated by the land side. This routine closes that gap: a lake-only geological unit
  !>                 that is missing from the geology LUT is reported here instead of failing deeper in the
  !>                 per-cell baseflow lookup.

  !     CALLING SEQUENCE
  !         call lake_check_geo_units(geo_unit, geo_unit_list)

  !     INTENT(IN)
  !>       \param[in] "integer(i4), dimension(:) :: geo_unit"       [1] geological unit ID on each footprint cell
  !>       \param[in] "integer(i4), dimension(:) :: geo_unit_list"  [1] geological unit IDs in the geology LUT

  !     HISTORY
  !>       \authors Matthias Zink, Rohini Kumar (land-surface geology-unit LUT check, from mo_mpr_container:mpr_check_geo_units_against_lut)
  !>       \date    January 2013
  !        Modifications:
  !        Pallav K Shrestha, October 2026 - applied the geology-unit presence check over the lake footprint
  !            (lake-only geological units are not validated by the land-surface check), for the detached lake MPR.

  subroutine lake_check_geo_units(geo_unit, geo_unit_list)

    implicit none

    integer(i4), dimension(:), intent(in) :: geo_unit       ! geological unit ID on each footprint cell
    integer(i4), dimension(:), intent(in) :: geo_unit_list  ! geological unit IDs in the geology LUT

    integer(i8) :: i

    do i = 1_i8, size(geo_unit, kind=i8)
      if (.not.any(geo_unit_list == geo_unit(i))) then
        log_fatal(*) "Lake MPR: lake footprint geological unit ", n2s(geo_unit(i)), &
          " is not present in the geology LUT."
        error stop 1
      end if
    end do

  end subroutine lake_check_geo_units


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
