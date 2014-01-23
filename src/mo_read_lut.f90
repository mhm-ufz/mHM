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

  PUBLIC :: read_gauge_lut         ! Reads LUT containing gauge information
  PUBLIC :: read_geoformation_lut  ! Reads LUT containing geological formation information
  PUBLIC :: read_lai_lut           ! Reads LUT containing LAI class information

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !         read_gauge_lut

  !     PURPOSE
  !>        \brief Reading gauge configurations for mHM.

  !>        \details This routine reads the setup which gauge is used for the
  !>                 calibration / evaluation of which catchment. 
  !>                 Please specify the ID of  the gauging stations within this file. 
  !>                 The ID has to correspond to the ID's given in the 'gaugelocation.asc' and
  !>                 to the filename containing the time series (the filename of the time series is
  !>                 created intrinsicly like "gaugeId.dat", e.g. gaugeID=000111 --> filename='000111.dat')
  !>                 please be aware that zeros are considered
  !> 
  !>                 structure:
  !>                 total number of gauges (all subasins, all basins) to be used for model evaluation
  !>                 header
  !>                 data corresponding to header

  !     INTENT(IN)
  !>        \param[in] "character(len=*)                :: filename"     Name of file
  !>        \param[in] "integer(i4),                    :: fileunit"     Unit of file to open
  !>        \param[in] "character(len=*),dimension(:),  :: dirGauges"    Directory of gauging data
  !>        \param[in] "integer(i4),                    :: nBasins"      Number of basins for simulation
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "integer(i4),    dimension(:),   allocatable :: gauge_basinId"        Basin Id
  !>        \param[out] "integer(i4),    dimension(:),   allocatable :: gauge_gaugeId"        Gauge-Id (e.g. 0000444)
  !>        \param[out] "character(256), dimension(:),   allocatable :: gauge_fname"          Name of runoff file
  !>        \param[out] "integer(i4),    dimension(:),   allocatable :: basin_nGauges"        Number of gauges within a basin
  !>        \param[out] "integer(i4),    dimension(:,:), allocatable :: basin_gaugeIdList"    Gauge Id list
  !>                                                                                          (e.g. 0000444 0000445)
  !>        \param[out] "integer(i4),    dimension(:,:), allocatable :: basin_gaugeIndexList" Gauge index list
  !>                                                                                          (e.g. 1 for 00444, 2 for 00445)
  !>        \param[out] "integer(i4),    dimension(:,:), allocatable :: basin_gaugeNodeList"  Node id of gauge -> netNode

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
  !>        \author Matthias Zink
  !>        \date Jan 2013

  ! ------------------------------------------------------------------

  SUBROUTINE read_gauge_lut(filename, fileunit, dirGauges, nBasins, &      ! Intent IN
       nGaugesTotal, gauge_basinId, gauge_gaugeId, gauge_fname,     &      ! Intent OUT
       basin_nGauges, basin_gaugeIdList, basin_gaugeIndexList,      &      ! Intent OUT
       basin_gaugeNodeList)                                                ! Intent OUT

    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_mhm_constants,    only: nodata_i4, maxNoGauges

    implicit none

    character(len=*),                              intent(in)  :: filename             ! file name
    integer(i4),                                   intent(in)  :: fileunit             ! unit to open file
    character(len=*), dimension(:),                intent(in)  :: dirGauges            ! directory of gauge data
    integer(i4),                                   intent(in)  :: nBasins              ! number of basins
    integer(i4),                                   intent(out) :: nGaugesTotal         ! number of evaluation/gauging stations
    integer(i4),      dimension(:),   allocatable, intent(out) :: gauge_basinId        ! Basin Id
    integer(i4),      dimension(:),   allocatable, intent(out) :: gauge_gaugeId        ! gauge-Id (e.g. 0000444)
    character(256),   dimension(:),   allocatable, intent(out) :: gauge_fname          ! name runoff file
    integer(i4),      dimension(:),   allocatable, intent(out) :: basin_nGauges        ! number of gauges within a basin
    integer(i4),      dimension(:,:), allocatable, intent(out) :: basin_gaugeIdList    ! gauge Id list
    !                                                                                  ! (e.g. 0000444 0000445)
    integer(i4),      dimension(:,:), allocatable, intent(out) :: basin_gaugeIndexList ! Gauge index list
    !                                                                                  ! (e.g. 1 for 00444, 2 for 00445)
    integer(i4),      dimension(:,:), allocatable, intent(out) :: basin_gaugeNodeList  ! gauged node -> netNode
    !
    ! local variables
    integer(i4)                            :: i, j 
    integer(i4)                            :: basinId
    integer(i4)                            :: nGauges
    integer(i4)                            :: gaugeCounter, maxGaugePerBasin
    integer(i4),    dimension(maxNoGauges) :: gaugeId
    character(256)                         :: dummy
    character(256), dimension(maxNoGauges) :: gaugefNamePart

    open(unit=fileunit, file=filename, action='read', status='old')
    read(fileunit, *) dummy, nGaugesTotal
    read(fileunit, *) dummy

    allocate(gauge_gaugeId(nGaugesTotal)) 
    allocate(gauge_basinId(nGaugesTotal))
    allocate(gauge_fName(nGaugesTotal))

    if (nGaugesTotal .GT. maxNoGauges) then
       call message('***ERROR: read_gauge_lut: Total number of gauges is restricted to', num2str(maxNoGauges))         
       stop
    end if

    ! determine maximum number of gauges per basin -> for allocating array
    maxGaugePerBasin = 0_i4
    gaugecounter = 0_i4
    do i = 1, nBasins
       ! Id's have to be read as characters to create filename
       read(fileunit, *) basinId, nGauges, (gaugefNamePart(j), j=1,nGauges)
       do j = 1, nGauges
          gaugeCounter = gaugeCounter + 1_i4
          ! write to global variable gauge
          dummy=trim(dirGauges(i)) // trim(gaugefNamePart(j)) // trim('.txt')
          gauge_fname(gaugeCounter) = dummy
       end do
       maxGaugePerBasin = max(maxGaugePerBasin, nGauges)
    end do

    if (nGaugesTotal .NE. gaugeCounter) then
       call message('***ERROR: read_gauge_lut: Total number of gauges is not consistent with sum of gauges for single basins!')
       call message('File:', trim(filename))
       stop
    end if

    allocate(basin_nGauges       (nBasins                  ))
    allocate(basin_gaugeIdList   (nBasins, maxGaugePerBasin))
    allocate(basin_gaugeIndexList(nBasins, maxGaugePerBasin))
    allocate(basin_gaugeNodeList (nBasins, maxGaugePerBasin))

    basin_nGauges     = nodata_i4
    basin_gaugeIdList = nodata_i4

    basin_gaugeNodeList =  nodata_i4

    rewind(fileunit)

    read(fileunit, *) dummy, nGaugesTotal
    read(fileunit, *) dummy

    gaugecounter = 0_i4
    do i = 1, nBasins
       ! Id's have to be read as characters to create filename
       read(fileunit, *) basinId, nGauges, (gaugeId(j), j=1,nGauges)
       !
       if (basinId .NE. i) then 
          call message('***ERROR: ID of single basins are not continuous!')
          call message('FILE:', trim(filename))
          stop
       end if
       !
       do j = 1, nGauges
          gaugeCounter = gaugeCounter + 1_i4
          ! write to global variable gauge
          gauge_basinId(gaugeCounter) = i
          gauge_gaugeId(gaugeCounter) = gaugeId(j)
          basin_gaugeIdList(i,j)      = gaugeId(j)
          basin_gaugeIndexList(i,j)   = gaugecounter
       end do
       basin_nGauges(i) = nGauges
    end do

    close(fileunit)

  END SUBROUTINE read_gauge_lut

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
