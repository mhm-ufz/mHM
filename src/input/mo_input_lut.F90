!> \file mo_input_lut.F90
!> \copydoc mo_input_lut

!> \brief Routines reading lookup tables (lut).
!> \details This module contains routines reading various lookup tables (lut).
!! Reads soil records, top-horizon texture, geological formations and LAI classes.
!! Returned data belongs to the caller; regionalization and soil-horizon construction
!! remain in the consuming component.
!> \authors Juliane Mai, Matthias Zink, Rohini Kumar, Pallav K Shrestha, Sebastian Mueller
!> \date Jan 2013
!> \version 0.1
!> \changelog
!! - 2026: shared input decoding for land and lake regionalization; retain donor formats and raw values.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_input
#include "logging.h"
MODULE mo_input_lut

  ! Written    Juliane Mai,    Jan 2013
  ! Modified   Matthias Zink,  Jan 2013 - add read_gauge_lut

  use mo_logging
  USE mo_kind, ONLY : i4, dp
  USE mo_os, ONLY: check_path_isfile
  use mo_string_utils, ONLY: num2str
  use mo_message, ONLY: error_message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_geoformation_lut  ! Reads LUT containing geological formation information
  PUBLIC :: read_lai_lut           ! Reads LUT containing LAI class information

  public :: read_soil_lut_records, read_soil_texture_lut

  !> \class soil_lut_record_t
  !> \brief One raw soil-class record, before component-specific transformations.
  type, public :: soil_lut_record_t
    integer(i4) :: id = 0_i4 !< Soil class identifier.
    integer(i4) :: horizon = 1_i4 !< Horizon identifier; one for texture-only tables.
    real(dp) :: upper = 0.0_dp !< Upper depth in mm (mode 0).
    real(dp) :: lower = 0.0_dp !< Lower depth in mm (mode 0).
    real(dp) :: clay = 0.0_dp !< Clay content in percent.
    real(dp) :: sand = 0.0_dp !< Sand content in percent.
    real(dp) :: density = 0.0_dp !< Mineral bulk density in g/cm3.
  end type soil_lut_record_t

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        read_geoformation_lut

  !    PURPOSE
  !>       \brief Reads LUT containing geological formation information.

  !>       \details The LUT needs to have the following header:
  !>       \verbatim
  !>       nGeo_Formations  < Number of lines containing data >
  !>       GeoParam(i)   ClassUnit     Karstic      Description
  !>       \endverbatim

  !>       The subsequent lines contains the geological formation information:
  !>       \verbatim
  !>       <GeoParam(i)>  <ClassUnit_i4>  <Karstic_i4>  <Description_char>
  !>       \endverbatim
  !>       All following lines will be discarded while reading.
  !>       GeoParam(i) is the global GeoParam table index while ClassUnit is the unit of the map containing the
  !>       geological formations. ClassUnit values do not necessarily contain subsequent numbers. The parametrization
  !>       of this unit is part of the `baseflow_1` parameter namelist.

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Jan 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_geoformation_lut(filename, fileunit, nGeo, geo_unit, geo_karstic, geo_param_index)
    implicit none
    character(len = *), intent(in) :: filename !< Input lookup-table path.
    integer(i4), intent(in), optional :: fileunit !< Optional unit to open file
    integer(i4), intent(out) :: nGeo !< Number of geological formations
    integer(i4), dimension(:), allocatable, intent(out) :: geo_unit !< List of id numbers of each geological formations
    integer(i4), dimension(:), allocatable, intent(out) :: geo_karstic !< ID of the Karstic formation (0 == does not exist)
    integer(i4), dimension(:), allocatable, intent(out), optional :: geo_param_index !< Global GeoParam index for each geological formation

    integer(i4) :: i, ios, unit
    integer(i4) :: param_index

    character(256) :: dummy

    if (present(fileunit)) then
      unit = fileunit
      open(unit=unit, file=filename, action='read', status='old', iostat=ios)
    else
      open(newunit=unit, file=filename, action='read', status='old', iostat=ios)
    end if
    if (ios /= 0) then
      log_fatal(*) "LUT: cannot open geology table: ", trim(filename)
      error stop 1
    end if
    read(unit, *, iostat=ios) dummy, nGeo
    if (ios /= 0) then
      log_fatal(*) "LUT: malformed geology table header: ", trim(filename)
      error stop 1
    end if
    if (nGeo < 1_i4) then
      log_fatal(*) "LUT: geology table must declare at least one class: ", trim(filename)
      error stop 1
    end if
    read(unit, *, iostat=ios) dummy
    if (ios /= 0) then
      log_fatal(*) "LUT: missing geology column header: ", trim(filename)
      error stop 1
    end if
    if (present(geo_param_index)) allocate(geo_param_index(nGeo))
    allocate(geo_unit(nGeo), geo_karstic(nGeo))

    ! read data
    do i = 1, nGeo
      read(unit, *, iostat=ios) param_index, geo_unit(i), geo_karstic(i), dummy
      if (ios /= 0) then
        log_fatal(*) "LUT: declared geology classes exceed available valid rows: ", trim(filename)
        error stop 1
      end if
      if (present(geo_param_index)) geo_param_index(i) = param_index
    end do

    close(unit)

  end subroutine read_geoformation_lut

  ! ------------------------------------------------------------------

  !    NAME
  !        read_lai_lut

  !    PURPOSE
  !>       \brief Reads LUT containing LAI information.

  !>       \details The LUT needs to have the following header:
  !>       \verbatim
  !>       NoLAIclasses  <Number of lines containing data>
  !>       Id  land-use  Jan.   Feb.    Mar.    Apr.    May    Jun.    Jul.    Aug.    Sep.    Oct.    Nov.    Dec.
  !>       \endverbatim
  !>       The subsequent lines contains the lai class information:
  !>       \verbatim
  !>       <ID_i4>  <landuse_char>  <val_1_dp>  <val_2_dp>  <val_3_dp>  <val_4_dp> ... <val_12_dp>
  !>       \endverbatim
  !>       All following lines will be discarded while reading.

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Jan 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_lai_lut(filename, fileunit, nLAI, LAIIDlist, LAI)

    use mo_constants, only : YearMonths

    implicit none
    character(len = *), intent(in) :: filename !< Input lookup-table path.
    integer(i4), intent(in), optional :: fileunit !< Optional unit to open file
    integer(i4), intent(out) :: nLAI !< Number of LAI classes
    integer(i4), dimension(:), allocatable, intent(out) :: LAIIDlist !< List of ids of LAI classes
    real(dp), dimension(:, :), allocatable, intent(out) :: LAI !< LAI per class (row) and month (col)

    integer(i4) :: i, j, ios, unit

    character(256) :: dummy

    !checking whether the file exists
    call check_path_isfile(path = filename, raise=.true.)
    if (present(fileunit)) then
      unit = fileunit
      open(unit=unit, file=filename, action='read', status='old')
    else
      open(newunit=unit, file=filename, action='read', status='old')
    end if

    ! read header
    read(unit, *) dummy, nLAI
    read(unit, *) dummy
    dummy = dummy // ''   ! only to avoid warning

    ! allocate arrays
    allocate(LAIIDList(nLAI))
    allocate(LAI(nLAI, int(YearMonths, i4)))

    ! read data
    do i = 1, nLAI
      read(unit, *, iostat=ios) LAIIDList(i), dummy, (LAI(i, j), j = 1, int(YearMonths, i4))
      if ( ios /= 0 ) call error_message( &
        "ERROR: NoLAIclasses (", num2str(nLAI), ") in LAI_classdefinition.txt ", &
        "seems to be higher than available LAI classes!" &
      )
    end do

    close(unit)

  end subroutine read_lai_lut

  !> \brief Decode either soil LUT layout without changing physical values.
  !> \details Format parsing derives from Juliane Mai's read_soil_LUT and the lake
  !! texture reader adapted by Pallav K Shrestha. Consumers own returned records.
  subroutine read_soil_lut_records(path, n_classes, mode, records)
    character(*), intent(in) :: path !< Soil LUT path.
    integer(i4), intent(out) :: n_classes !< Declared soil-class count.
    integer(i4), intent(out) :: mode !< Detected layout: zero for horizons, one for textures.
    type(soil_lut_record_t), allocatable, intent(out) :: records(:) !< Raw rows in file order.
    type(soil_lut_record_t) :: record
    integer :: unit, ios, row, count
    character(1024) :: line
    character(256) :: header

    open(newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      log_fatal(*) "LUT: cannot open soil table: ", trim(path)
      error stop 1
    end if
    read(unit, *, iostat=ios) header, n_classes
    if (ios /= 0) then
      log_fatal(*) "LUT: malformed soil table header: ", trim(path)
      error stop 1
    end if
    if (n_classes < 1_i4) then
      log_fatal(*) "LUT: soil table must declare at least one class: ", trim(path)
      error stop 1
    end if
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) then
      log_fatal(*) "LUT: missing soil column header: ", trim(path)
      error stop 1
    end if
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios > 0) then
        log_fatal(*) "LUT: cannot read soil table: ", trim(path)
        error stop 1
      end if
      if (len_trim(line) > 0) count = count + 1
    end do
    if (count == 0) then
      log_fatal(*) "LUT: soil table has no records: ", trim(path)
      error stop 1
    end if
    allocate(records(count))
    rewind(unit)
    read(unit, *)
    read(unit, *)
    mode = -1_i4
    row = 0
    do while (row < count)
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) then
        log_fatal(*) "LUT: cannot reread soil table: ", trim(path)
        error stop 1
      end if
      if (len_trim(line) == 0) cycle
      if (mode < 0_i4) then
        read(line, *, iostat=ios) record%id, record%horizon, record%upper, record%lower, &
          record%clay, record%sand, record%density
        mode = 1_i4
        if (ios == 0) mode = 0_i4
      end if
      if (mode == 0_i4) then
        read(line, *, iostat=ios) record%id, record%horizon, record%upper, record%lower, &
          record%clay, record%sand, record%density
      else
        record%horizon = 1_i4
        read(line, *, iostat=ios) record%id, record%clay, record%sand, record%density
      end if
      if (ios /= 0) then
        log_fatal(*) "LUT: malformed soil table row: ", trim(path)
        error stop 1
      end if
      if (record%id < 1_i4 .or. record%id > n_classes .or. record%horizon < 1_i4) then
        log_fatal(*) "LUT: soil class or horizon outside declared range: ", trim(path)
        error stop 1
      end if
      row = row + 1
      records(row) = record
    end do
    close(unit)
  end subroutine read_soil_lut_records

  !     NAME
  !         read_soil_texture_lut

  !     PURPOSE
  !>        \brief Reads the top-horizon sand and clay content per soil class from the soil LUT.

  !>        \details Parses the mHM soil class-definition look-up table and returns, for every soil
  !>                 class, the sand and clay content of its uppermost horizon. Both soil_db_mode
  !>                 layouts of mo_soil_database:read_soil_LUT are supported, detected from the first
  !>                 data row and latched: a first line "nSoil_Types <n>", a header line, then one row
  !>                 per class in either mode 0 (class/horizon, 7 columns) "<class> <horizon>
  !>                 <upper_depth> <lower_depth> <clay%> <sand%> <bulk_density>" or mode 1 (texture,
  !>                 4 columns) "<class> <clay%> <sand%> <bulk_density>". In mode 0 only horizon 1 is
  !>                 retained; mode 1 carries one texture per class (its single horizon). Classes
  !>                 absent from the table keep the missing-value sentinel. Raw values are shared with the land reader;
  !>                 component-specific soil-database transformations are not applied here.

  !     CALLING SEQUENCE
  !         call read_soil_texture_lut(path, n_classes, sand_by_class, clay_by_class)

  !     HISTORY
  !>       \authors Juliane Mai (soil look-up table parsing, from mo_soil_database:read_soil_LUT)
  !>       \date    December 2012
  !        Modifications:
  !        Pallav K Shrestha, October 2026 - reused the soil-LUT parse logic as a self-contained reader
  !            returning only the top-horizon sand/clay per class, for the detached lake MPR.

  subroutine read_soil_texture_lut(path, n_classes, sand_by_class, clay_by_class)

    use mo_constants, only : nodata_dp

    implicit none

    character(*),              intent(in) :: path !< soil LUT path
    integer(i4),               intent(out) :: n_classes !< number of soil classes
    real(dp), allocatable,     intent(out) :: sand_by_class(:) !< [%] horizon-1 sand per class
    real(dp), allocatable,     intent(out) :: clay_by_class(:) !< [%] horizon-1 clay per class

    type(soil_lut_record_t), allocatable :: records(:)
    integer(i4) :: mode, row, id

    call read_soil_lut_records(path, n_classes, mode, records)
    allocate(sand_by_class(n_classes), source=nodata_dp)
    allocate(clay_by_class(n_classes), source=nodata_dp)
    do row = 1_i4, size(records)
      if (records(row)%horizon /= 1_i4) cycle
      id = records(row)%id
      sand_by_class(id) = records(row)%sand
      clay_by_class(id) = records(row)%clay
    end do
  end subroutine read_soil_texture_lut

END MODULE mo_input_lut
