!>       \file mo_read_lut.f90

!>       \brief Routines reading lookup tables (lut).

!>       \details This module contains routines reading various lookup tables (lut).
!>       (1) LUT containing gauge information.
!>       (2) LUT containing geological formation information.
!>       (3) LUT containing LAI class information.

!>       \authors Juliane Mai, Matthias Zink

!>       \date Jan 2013

! Modifications:

MODULE mo_read_lut

  ! Written    Juliane Mai,    Jan 2013
  ! Modified   Matthias Zink,  Jan 2013 - add read_gauge_lut

  USE mo_kind, ONLY : i4, dp
  USE mo_os, ONLY: path_isfile

  IMPLICIT NONE

  PUBLIC :: read_geoformation_lut  ! Reads LUT containing geological formation information
  PUBLIC :: read_lai_lut           ! Reads LUT containing LAI class information

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
  !>       GeoParam is a running index while ClassUnit is the unit of the map containing the geological formations
  !>       such that it does not neccessarily contains subsequent numbers. The parametrization of this unit is part
  !>       of the namelist mhm_parameter.nml under <geoparameter>.

  !    INTENT(IN)
  !>       \param[in] "character(len = *) :: filename" File name of LUT
  !>       \param[in] "integer(i4) :: fileunit"        Unit to open file

  !    INTENT(OUT)
  !>       \param[out] "integer(i4) :: nGeo"                      Number of geological formations
  !>       \param[out] "integer(i4), dimension(:) :: geo_unit"    List of id numbers of each geological formations
  !>       \param[out] "integer(i4), dimension(:) :: geo_karstic" ID of the Karstic formation (0 == does not exist)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Jan 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_geoformation_lut(filename, fileunit, nGeo, geo_unit, geo_karstic)
    implicit none

    ! File name of LUT
    character(len = *), intent(in) :: filename

    ! Unit to open file
    integer(i4), intent(in) :: fileunit

    ! Number of geological formations
    integer(i4), intent(out) :: nGeo

    ! List of id numbers of each geological formations
    integer(i4), dimension(:), allocatable, intent(out) :: geo_unit

    ! ID of the Karstic formation (0 == does not exist)
    integer(i4), dimension(:), allocatable, intent(out) :: geo_karstic

    integer(i4) :: i

    character(256) :: dummy

    !checking whether the file exists
    call path_isfile(path = filename, quiet_ = .true., throwError_ = .true.)
    open(fileunit, file = filename, action = 'read', status = 'old')

    ! read header
    read(fileunit, *) dummy, nGeo
    read(fileunit, *) dummy
    dummy = dummy // ''   ! only to avoid warning

    ! allocation of arrays
    allocate(geo_unit(nGeo))
    allocate(geo_karstic(nGeo))

    ! read data
    do i = 1, nGeo
      read(fileunit, *) dummy, geo_unit(i), geo_karstic(i), dummy
    end do

    close(fileunit)

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

  !    INTENT(IN)
  !>       \param[in] "character(len = *) :: filename" File name of LUT
  !>       \param[in] "integer(i4) :: fileunit"        Unit to open file

  !    INTENT(OUT)
  !>       \param[out] "integer(i4) :: nLAIs"                    Number of LAI classes
  !>       \param[out] "integer(i4), dimension(:) :: LAIIDlist" List of ids of LAI classes
  !>       \param[out] "real(dp), dimension(:, :) :: LAI"       LAI per class (row) and month (col)

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Jan 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine read_lai_lut(filename, fileunit, nLAIs, LAIIDlist, LAI)

    use mo_constants, only : YearMonths

    implicit none

    ! File name of LUT
    character(len = *), intent(in) :: filename

    ! Unit to open file
    integer(i4), intent(in) :: fileunit

    ! Number of LAI classes
    integer(i4), intent(out) :: nLAIs

    ! List of ids of LAI classes
    integer(i4), dimension(:), allocatable, intent(out) :: LAIIDlist

    ! LAI per class (row) and month (col)
    real(dp), dimension(:, :), allocatable, intent(out) :: LAI

    integer(i4) :: i, j

    character(256) :: dummy

    !checking whether the file exists
    call path_isfile(path = filename, quiet_ = .true., throwError_ = .true.)
    open(fileunit, file = filename, action = 'read')

    ! read header
    read(fileunit, *) dummy, nLAIs
    read(fileunit, *) dummy
    dummy = dummy // ''   ! only to avoid warning

    ! allocate arrays
    allocate(LAIIDList(nLAIs))
    allocate(LAI(nLAIs, nint(YearMonths, i4)))

    ! read data
    do i = 1, nLAIs
      read(fileunit, *) LAIIDList(i), dummy, (LAI(i, j), j = 1, nint(YearMonths, i4))
    end do

    close(fileunit)

  end subroutine read_lai_lut

END MODULE mo_read_lut
