!> \file mo_read_lut.f90

!> \brief Routines reading lookup tables (lut).

!> \details This module contains routines reading various lookup tables (lut).\n
!> (1) LUT containing gauge information.\n
!> (2) LUT containing geological formation information.\n
!> (3) LUT containing LAI class information.\n

!> \authors Juliane Mai, Matthias Zink
!> \date Jan 2013

MODULE mo_read_lut

  ! Written    Juliane Mai,    Jan 2013
  ! Modified   Matthias Zink,  Jan 2013 - add read_gauge_lut

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PUBLIC :: read_geoformation_lut  ! Reads LUT containing geological formation information
  PUBLIC :: read_lai_lut           ! Reads LUT containing LAI class information

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !         read_geoformation_lut

  !     PURPOSE
  !>        \brief Reads LUT containing geological formation information.

  !>        \details The LUT needs to have the following header:\n
  !>        \verbatim
  !>        nGeo_Formations  < Number of lines containing data > 
  !>        GeoParam(i)   ClassUnit     Karstic      Description
  !>        \endverbatim
  !>        
  !>        The subsequent lines contains the geological formation information:\n
  !>        \verbatim
  !>        <GeoParam(i)>  <ClassUnit_i4>  <Karstic_i4>  <Description_char>
  !>        \endverbatim
  !>        All following lines will be discarded while reading.\n
  !>        GeoParam is a running index while ClassUnit is the unit of the map containing the geological formations 
  !>        such that it does not neccessarily contains subsequent numbers. The parametrization of this unit is part 
  !>        of the namelist mhm_parameter.nml under <geoparameter>.


  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: filename"        File name of LUT
  !>        \param[in] "integer(i4)      :: fileunit"        Unit to open file

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4)                              :: nGeo"           
  !>                                            Number of geological formations
  !>        \param[out] "integer(i4), dimension(:),   allocatable :: geo_unit"
  !>                                            List of id numbers of each geological formations
  !>        \param[out] "integer(i4), dimension(:),   allocatable :: geo_karstic"
  !>                                            ID of the Karstic formation (0 == does not exist)

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jan 2013

  ! ------------------------------------------------------------------

  subroutine read_geoformation_lut(filename, fileunit, nGeo, geo_unit, geo_karstic)

    implicit none 

    character(len=*),                         intent(in)  :: filename      ! name of file
    integer(i4),                              intent(in)  :: fileunit      ! unit to open file
    integer(i4),                              intent(out) :: nGeo          ! number of geological formations
    integer(i4), dimension(:),   allocatable, intent(out) :: geo_unit      ! list of id numbers of each geological formations
    integer(i4), dimension(:),   allocatable, intent(out) :: geo_karstic   ! id of the Karstic formation (0 == does not exist)

    ! local variables
    integer(i4)    :: i
    character(256) :: dummy

    open(fileunit, file=filename, action='read', status='old')

    ! read header
    read(fileunit, *) dummy, nGeo
    read(fileunit, *) dummy
    dummy = dummy//''   ! only to avoid warning

    ! allocation of arrays
    allocate(geo_unit(nGeo))
    allocate(geo_karstic(nGeo))

    ! read data
    do i=1,nGeo
       read(fileunit,*) dummy, geo_unit(i), geo_karstic(i), dummy
    end do

    close(fileunit)

  end subroutine read_geoformation_lut

  ! ------------------------------------------------------------------

  !      NAME
  !         read_lai_lut

  !     PURPOSE
  !>        \brief Reads LUT containing LAI information.

  !>        \details The LUT needs to have the following header:\n
  !>        \verbatim
  !>        NoLAIclasses  <Number of lines containing data>
  !>        Id  land-use  Jan.   Feb.    Mar.    Apr.    May    Jun.    Jul.    Aug.    Sep.    Oct.    Nov.    Dec.
  !>        \endverbatim
  !>        The subsequent lines contains the lai class information:\n
  !>        \verbatim
  !>        <ID_i4>  <landuse_char>  <val_1_dp>  <val_2_dp>  <val_3_dp>  <val_4_dp> ... <val_12_dp>
  !>        \endverbatim
  !>        All following lines will be discarded while reading.

  !     INTENT(IN)
  !>        \param[in] "character(len=*) :: filename"        File name of LUT
  !>        \param[in] "integer(i4)      :: fileunit"        Unit to open file

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4)                              :: nLAI"      Number of LAI classes
  !>        \param[out] "integer(i4), dimension(:),   allocatable :: LAIIDlist" List of ids of LAI classes
  !>        \param[out] "real(dp),    dimension(:,:), allocatable :: LAI"       LAI per class (row) and month (col)

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Jan 2013

  ! ------------------------------------------------------------------

  subroutine read_lai_lut(filename, fileunit, nLAI, LAIIDlist, LAI)

    use mo_mhm_constants, only: YearMonths  ! months per year

    implicit none 

    character(len=*),                         intent(in)  :: filename      ! name of file
    integer(i4),                              intent(in)  :: fileunit      ! unit to open file
    integer(i4),                              intent(out) :: nLAI          ! number of LAI classes
    integer(i4), dimension(:),   allocatable, intent(out) :: LAIIDlist     ! List of ids of LAI classes
    real(dp),    dimension(:,:), allocatable, intent(out) :: LAI           ! LAI per class (row) and month (col)

    ! local variables
    integer(i4)    :: i, j
    character(256) :: dummy

    open(fileunit, file=filename, action='read')

    ! read header
    read(fileunit, *) dummy, nLAI
    read(fileunit, *) dummy
    dummy = dummy//''   ! only to avoid warning

    ! allocate arrays
    allocate(LAIIDList(nLAI))
    allocate(LAI(nLAI, int(YearMonths,i4)))

    ! read data
    do i=1,nLAI
       read(fileunit, *) LAIIDList(i), dummy, (LAI(i,j),j=1,int(YearMonths,i4))
    end do

    close(fileunit)

  end subroutine read_lai_lut

END MODULE mo_read_lut
