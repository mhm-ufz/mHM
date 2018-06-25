!>       \file mo_mrm_read_config.f90

!>       \brief read mRM config

!>       \details This module contains all mRM subroutines related to
!>       reading the mRM configuration either from file or copy from mHM.

!>       \authors Stephan Thober

!>       \date Aug 2015

! Modifications:

module mo_mrm_read_config

  use mo_kind, only : i4, dp

  implicit none

  public :: mrm_read_config

contains
  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_read_config

  !    PURPOSE
  !>       \brief Read the general config of mRM

  !>       \details Depending on the variable mrm_coupling_config, the
  !>       mRM config is either read from mrm.nml and parameters from
  !>       mrm_parameter.nml or copied from mHM.

  !    INTENT(IN)
  !>       \param[in] "character(*) :: file_namelist, file_namelist_param"
  !>       \param[in] "integer :: unamelist, unamelist_param"
  !>       \param[in] "character(*) :: file_namelist, file_namelist_param"
  !>       \param[in] "integer :: unamelist, unamelist_param"
  !>       \param[in] "logical :: do_message"                              - flag for writing mHM standard messages

  !    INTENT(OUT)
  !>       \param[out] "logical :: readLatLon" - flag for reading LatLon file

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober  Sep 2015 - removed stop condition when routing resolution is smaller than hydrologic resolution
  ! Stephan Thober  Oct 2015 - added NLoutputResults namelist, fileLatLon to directories_general namelist, and readLatLon flag
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_config(file_namelist, unamelist, file_namelist_param, unamelist_param, do_message, readLatLon)

    use mo_common_constants, only : maxNoBasins, nodata_i4
    use mo_common_mHM_mRM_read_config, only : common_check_resolution
    use mo_common_variables, only : ALMA_convention, nBasins, processMatrix
    use mo_message, only : message
    use mo_mrm_constants, only : maxNoGauges
    use mo_mrm_file, only : file_defOutput, udefOutput
    use mo_mrm_global_variables, only : InflowGauge, basinInfo_mRM, basin_mrm, &
                                        dirGauges, dirTotalRunoff, filenameTotalRunoff, gauge, is_start, nGaugesTotal, &
                                        nInflowGaugesTotal, outputFlxState_mrm, timeStep_model_outputs_mrm, &
                                        varnameTotalRunoff
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str

    implicit none

    character(*), intent(in) :: file_namelist, file_namelist_param

    integer, intent(in) :: unamelist, unamelist_param

    ! - flag for writing mHM standard messages
    logical, intent(in) :: do_message

    ! - flag for reading LatLon file
    logical, intent(out) :: readLatLon

    integer(i4), dimension(maxNoBasins) :: NoGauges_basin

    integer(i4), dimension(maxNoBasins, maxNoGauges) :: Gauge_id

    character(256), dimension(maxNoBasins, maxNoGauges) :: Gauge_filename

    integer(i4), dimension(maxNoBasins) :: NoInflowGauges_basin

    integer(i4), dimension(maxNoBasins, maxNoGauges) :: InflowGauge_id

    character(256), dimension(maxNoBasins, maxNoGauges) :: InflowGauge_filename

    logical, dimension(maxNoBasins, maxNoGauges) :: InflowGauge_Headwater

    integer(i4) :: iBasin

    integer(i4) :: iGauge

    integer(i4) :: idx

    character(256), dimension(maxNoBasins) :: dir_Gauges

    character(256), dimension(maxNoBasins) :: dir_Total_Runoff

    logical :: file_exists

    type(basinInfo_mRM), pointer :: basin_mrm_iBasin


    ! namelist spatial & temporal resolution, optmization information
    namelist /mainconfig_mrm/ ALMA_convention, filenameTotalRunoff, varnameTotalRunoff
    ! namelist directories
    namelist /directories_mRM/ dir_Gauges, dir_Total_Runoff
    namelist /evaluation_gauges/ nGaugesTotal, NoGauges_basin, Gauge_id, gauge_filename
    ! namelist for inflow gauges
    namelist /inflow_gauges/ nInflowGaugesTotal, NoInflowGauges_basin, InflowGauge_id, &
            InflowGauge_filename, InflowGauge_Headwater
    ! name list regarding output
    namelist /NLoutputResults/timeStep_model_outputs_mrm, outputFlxState_mrm

    !===============================================================
    ! INITIALIZATION
    !===============================================================
    is_start = .True.
    nGaugesTotal = nodata_i4
    NoGauges_basin = nodata_i4
    Gauge_id = nodata_i4
    gauge_filename = num2str(nodata_i4)

    ! default arguments
    ALMA_convention = .false.
    filenameTotalRunoff = 'total_runoff'
    varnameTotalRunoff = 'total_runoff'

    !===============================================================
    !  Read namelist main directories
    !===============================================================
    call open_nml(file_namelist, unamelist, quiet = .true.)

    !===============================================================
    !  Read namelist for mainconfig for mRM
    !===============================================================
    call position_nml('mainconfig_mrm', unamelist)
    read(unamelist, nml = mainconfig_mrm)

    !===============================================================
    !  Read namelist for mainpaths
    !===============================================================
    call position_nml('directories_mRM', unamelist)
    read(unamelist, nml = directories_mRM)

    dirGauges = dir_Gauges(1 : nBasins)
    dirTotalRunoff = dir_Total_Runoff(1 : nBasins)

    !===============================================================
    ! READ EVALUATION GAUGES
    !===============================================================
    call position_nml('evaluation_gauges', unamelist)
    read(unamelist, nml = evaluation_gauges)

    if (nGaugesTotal .GT. maxNoGauges) then
      call message()
      call message('***ERROR: ', trim(file_namelist), ': Total number of evaluation gauges is restricted to', &
              num2str(maxNoGauges))
      call message('          Error occured in namlist: evaluation_gauges')
      stop 1
    end if

    allocate(gauge%gaugeId(nGaugesTotal)) ; gauge%gaugeId = nodata_i4
    allocate(gauge%basinId(nGaugesTotal)) ; gauge%basinId = nodata_i4
    allocate(gauge%fName  (nGaugesTotal)) ; gauge%fName(1) = num2str(nodata_i4)
    allocate(basin_mrm(nBasins))

    idx = 0
    do iBasin = 1, nBasins
      basin_mrm_iBasin => basin_mrm(iBasin)
      ! initialize
      basin_mrm_iBasin%nGauges = nodata_i4
      allocate(basin_mrm_iBasin%gaugeIdList(maxval(NoGauges_basin(:))))
      basin_mrm_iBasin%gaugeIdList = nodata_i4
      allocate(basin_mrm_iBasin%gaugeIndexList(maxval(NoGauges_basin(:))))
      basin_mrm_iBasin%gaugeIndexList = nodata_i4
      allocate(basin_mrm_iBasin%gaugeNodeList(maxval(NoGauges_basin(:))))
      basin_mrm_iBasin%gaugeNodeList = nodata_i4
      ! check if NoGauges_basin has a valid value
      if (NoGauges_basin(iBasin) .EQ. nodata_i4) then
        call message()
        call message('***ERROR: ', trim(file_namelist), ': Number of evaluation gauges for subbasin ', &
                trim(adjustl(num2str(iBasin))), ' is not defined!')
        call message('          Error occured in namelist: evaluation_gauges')
        stop 1
      end if

      basin_mrm_iBasin%nGauges = NoGauges_basin(iBasin)

      do iGauge = 1, NoGauges_basin(iBasin)
        ! check if NoGauges_basin has a valid value
        if (Gauge_id(iBasin, iGauge) .EQ. nodata_i4) then
          call message()
          call message('***ERROR: ', trim(file_namelist), ': ID ', &
                  trim(adjustl(num2str(Gauge_id(iBasin, iGauge)))), ' of evaluation gauge ', &
                  trim(adjustl(num2str(iGauge))), ' for subbasin ', &
                  trim(adjustl(num2str(iBasin))), ' is not defined!')
          call message('          Error occured in namelist: evaluation_gauges')
          stop 1
        else if (trim(gauge_filename(iBasin, iGauge)) .EQ. trim(num2str(nodata_i4))) then
          call message()
          call message('***ERROR: ', trim(file_namelist), ': Filename of evaluation gauge ', &
                  trim(adjustl(num2str(iGauge))), ' for subbasin ', &
                  trim(adjustl(num2str(iBasin))), ' is not defined!')
          call message('          Error occured in namelist: evaluation_gauges')
          stop 1
        end if
        !
        idx = idx + 1
        gauge%basinId(idx) = iBasin
        gauge%gaugeId(idx) = Gauge_id(iBasin, iGauge)
        gauge%fname(idx) = trim(dirGauges(iBasin)) // trim(gauge_filename(iBasin, iGauge))
        basin_mrm_iBasin%gaugeIdList(iGauge) = Gauge_id(iBasin, iGauge)
        basin_mrm_iBasin%gaugeIndexList(iGauge) = idx
      end do
    end do

    if (nGaugesTotal .NE. idx) then
      call message()
      call message('***ERROR: ', trim(file_namelist), ': Total number of evaluation gauges (', &
              trim(adjustl(num2str(nGaugesTotal))), &
              ') different from sum of gauges in subbasins (', trim(adjustl(num2str(idx))), ')!')
      call message('          Error occured in namelist: evaluation_gauges')
      stop
    end if

    !===============================================================
    ! Read inflow gauge information
    !===============================================================

    nInflowGaugesTotal = 0
    NoInflowGauges_basin = 0
    InflowGauge_id = nodata_i4
    InflowGauge_filename = num2str(nodata_i4)

    call position_nml('inflow_gauges', unamelist)
    read(unamelist, nml = inflow_gauges)

    if (nInflowGaugesTotal .GT. maxNoGauges) then
      call message()
      call message('***ERROR: ', trim(file_namelist), &
              ':read_gauge_lut: Total number of inflow gauges is restricted to', num2str(maxNoGauges))
      call message('          Error occured in namlist: inflow_gauges')
      stop
    end if

    ! allocation - max() to avoid allocation with zero, needed for mhm call
    allocate(InflowGauge%gaugeId (max(1, nInflowGaugesTotal)))
    allocate(InflowGauge%basinId (max(1, nInflowGaugesTotal)))
    allocate(InflowGauge%fName   (max(1, nInflowGaugesTotal)))
    ! dummy initialization
    InflowGauge%gaugeId = nodata_i4
    InflowGauge%basinId = nodata_i4
    InflowGauge%fName = num2str(nodata_i4)

    idx = 0
    do iBasin = 1, nBasins
      basin_mrm_iBasin => basin_mrm(iBasin)

      allocate(basin_mrm_iBasin%InflowGaugeIdList    (max(1, maxval(NoInflowGauges_basin(:)))))
      allocate(basin_mrm_iBasin%InflowGaugeHeadwater (max(1, maxval(NoInflowGauges_basin(:)))))
      allocate(basin_mrm_iBasin%InflowGaugeIndexList (max(1, maxval(NoInflowGauges_basin(:)))))
      allocate(basin_mrm_iBasin%InflowGaugeNodeList  (max(1, maxval(NoInflowGauges_basin(:)))))
      ! dummy initialization
      basin_mrm_iBasin%nInflowGauges = 0
      basin_mrm_iBasin%InflowGaugeIdList = nodata_i4
      basin_mrm_iBasin%InflowGaugeHeadwater = .FALSE.
      basin_mrm_iBasin%InflowGaugeIndexList = nodata_i4
      basin_mrm_iBasin%InflowGaugeNodeList = nodata_i4
      ! no inflow gauge for subbasin i
      if (NoInflowGauges_basin(iBasin) .EQ. nodata_i4) then
        NoInflowGauges_basin(iBasin) = 0
      end if

      basin_mrm_iBasin%nInflowGauges = NoInflowGauges_basin(iBasin)

      do iGauge = 1, NoInflowGauges_basin(iBasin)
        ! check if NoInflowGauges_basin has a valid value
        if (InflowGauge_id(iBasin, iGauge) .EQ. nodata_i4) then
          call message()
          call message('***ERROR: ', trim(file_namelist), ':ID of inflow gauge ', &
                  trim(adjustl(num2str(iGauge))), ' for subbasin ', &
                  trim(adjustl(num2str(iBasin))), ' is not defined!')
          call message('          Error occured in namlist: inflow_gauges')
          stop
        else if (trim(InflowGauge_filename(iBasin, iGauge)) .EQ. trim(num2str(nodata_i4))) then
          call message()
          call message('***ERROR: ', trim(file_namelist), ':Filename of inflow gauge ', &
                  trim(adjustl(num2str(iGauge))), ' for subbasin ', &
                  trim(adjustl(num2str(iBasin))), ' is not defined!')
          call message('          Error occured in namlist: inflow_gauges')
          stop
        end if
        !
        idx = idx + 1
        InflowGauge%basinId(idx) = iBasin
        InflowGauge%gaugeId(idx) = InflowGauge_id(iBasin, iGauge)
        InflowGauge%fname(idx) = trim(dirGauges(iBasin)) // trim(InflowGauge_filename(iBasin, iGauge))
        basin_mrm_iBasin%InflowGaugeIdList(iGauge) = InflowGauge_id(iBasin, iGauge)
        basin_mrm_iBasin%InflowGaugeHeadwater(iGauge) = InflowGauge_Headwater(iBasin, iGauge)
        basin_mrm_iBasin%InflowGaugeIndexList(iGauge) = idx
      end do
    end do

    if (nInflowGaugesTotal .NE. idx) then
      call message()
      call message('***ERROR: ', trim(file_namelist), ': Total number of inflow gauges (', &
              trim(adjustl(num2str(nInflowGaugesTotal))), &
              ') different from sum of inflow gauges in subbasins (', trim(adjustl(num2str(idx))), ')!')
      call message('          Error occured in namlist: inflow_gauges')
      stop
    end if

    call common_check_resolution(do_message, .true.)

    call close_nml(unamelist)

    !===============================================================
    ! Read namelist global parameters
    !===============================================================
    call read_mrm_routing_params(processMatrix(8, 1), file_namelist_param, unamelist_param)

    !===============================================================
    ! Read Output specifications for mRM
    !===============================================================
    outputFlxState_mrm = .FALSE.
    timeStep_model_outputs_mrm = -2
    inquire(file = file_defOutput, exist = file_exists)
    if (file_exists) then
      ! file exists
      call open_nml(file_defOutput, udefOutput, quiet = .true.)
      call position_nml('NLoutputResults', udefOutput)
      read(udefOutput, nml = NLoutputResults)
      call close_nml(udefOutput)
    else
      call message('')
      call message('No file specifying mRM output fluxes exists')
    end if
    readLatLon = any(outputFlxState_mrm)

    if (any(outputFlxState_mrm)) then
      call message('')
      call message('    Following output will be written:')

      call message('    FLUXES:')
      if (outputFlxState_mrm(1)) then
        call message('      routed streamflow      (L11_qMod)                [mm]')
      end if
    end if

    call message('')
    call message('    FINISHED reading config')
    call message('')

  end subroutine mrm_read_config

  ! ---------------------------------------------------------------------------
  ! SUBROUTINE READ_MRM_ROUTING_PARAMS
  ! ---------------------------------------------------------------------------
  !    NAME
  !        read_mrm_routing_params

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: processCase"          it is the default case, should be one
  !>       \param[in] "character(*) :: file_namelist_param" file name containing parameter namelist
  !>       \param[in] "integer(i4) :: unamelist_param"      file name id containing parameter namelist

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine read_mrm_routing_params(processCase, file_namelist_param, unamelist_param)

    use mo_common_constants, only : nColPars
    use mo_common_functions, only : in_bound
    use mo_common_variables, only : global_parameters, global_parameters_name, processMatrix
    use mo_message, only : message
    use mo_nml, only : close_nml, open_nml, position_nml
#ifndef MRM2MHM
    use mo_append, only : append
#endif

    implicit none

    ! it is the default case, should be one
    integer(i4), intent(in) :: processCase

    ! file name containing parameter namelist
    character(*), intent(in) :: file_namelist_param

    ! file name id containing parameter namelist
    integer(i4), intent(in) :: unamelist_param

#ifdef MRM2MHM
    ! equals sum of previous parameters
    integer(i4) :: start_index

#endif
    real(dp), dimension(nColPars) :: muskingumTravelTime_constant

    real(dp), dimension(nColPars) :: muskingumTravelTime_riverLength

    real(dp), dimension(nColPars) :: muskingumTravelTime_riverSlope

    real(dp), dimension(nColPars) :: muskingumTravelTime_impervious

    real(dp), dimension(nColPars) :: muskingumAttenuation_riverSlope

    real(dp), dimension(nColPars) :: streamflow_celerity
    real(dp), dimension(nColPars) :: g1
    real(dp), dimension(nColPars) :: g2

    namelist /routing1/ muskingumTravelTime_constant, muskingumTravelTime_riverLength, &
            muskingumTravelTime_riverSlope, muskingumTravelTime_impervious, muskingumAttenuation_riverSlope
    namelist /routing2/ streamflow_celerity
    namelist /routing3/ g1, g2
    !
    call open_nml(file_namelist_param, unamelist_param, quiet = .true.)

    if (processCase .eq. 1_i4) then
      call position_nml('routing1', unamelist_param)
      read(unamelist_param, nml = routing1)
    else if (processCase .eq. 2_i4) then
       call position_nml('routing2', unamelist_param)
       read(unamelist_param, nml = routing2)
    else if (processCase .eq. 3_i4) then
       call position_nml('routing3', unamelist_param)
       read(unamelist_param, nml = routing3)
    end if

#ifdef MRM2MHM
    ! -------------------------------------------------------------------------
    ! INCLUDE MRM PARAMETERS IN PARAMETERS OF MHM
    ! -------------------------------------------------------------------------
    ! Muskingum routing parameters with MPR
    if (processCase .eq. 1_i4) then
      ! insert parameter values and names at position required by mhm
      processMatrix(8, 1) = processCase
      processMatrix(8, 2) = 5_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
      start_index = processMatrix(8, 3) - processMatrix(8, 2)
      global_parameters(start_index + 1, :) = muskingumTravelTime_constant
      global_parameters(start_index + 2, :) = muskingumTravelTime_riverLength
      global_parameters(start_index + 3, :) = muskingumTravelTime_riverSlope
      global_parameters(start_index + 4, :) = muskingumTravelTime_impervious
      global_parameters(start_index + 5, :) = muskingumAttenuation_riverSlope

      global_parameters_name(start_index + 1 : start_index + processMatrix(8, 2)) = (/ &
              'muskingumTravelTime_constant   ', &
                      'muskingumTravelTime_riverLength', &
                      'muskingumTravelTime_riverSlope ', &
                      'muskingumTravelTime_impervious ', &
                      'muskingumAttenuation_riverSlope'/)
      ! adaptive timestep routing
    else if (processCase .eq. 2_i4) then
      processMatrix(8, 1) = processCase
      processMatrix(8, 2) = 1_i4
      processMatrix(8, 3) = sum(processMatrix(1 : 8, 2))
      start_index = processMatrix(8, 3) - processMatrix(8, 2)
      global_parameters(start_index + 1, :) = streamflow_celerity

      global_parameters_name(start_index + 1 : start_index + processMatrix(8, 2)) = (/ &
              'streamflow_celerity'/)
     ! adaptive timestep routing - varying celerity
     else if (processCase .eq. 3_i4) then
        ! insert parameter values and names at position required by mhm
        processMatrix(8, 1) = processCase
        processMatrix(8, 2) = 2_i4
        processMatrix(8, 3) = sum(processMatrix(1:8, 2))
        start_index         = processMatrix(8, 3) - processMatrix(8, 2)
        global_parameters(start_index + 1, :) = g1
        global_parameters(start_index + 2, :) = g2
    
        global_parameters_name(start_index + 1 : start_index + processMatrix(8,2)) = (/ &
             'g1', &
             'g2'/)
    end if
#else
    ! Muskingum routing parameters with MPR
    if (processCase .eq. 1_i4) then
      processMatrix(8, 1) = processCase
      processMatrix(8, 2) = 5_i4
      processMatrix(8, 3) = processMatrix(8, 2)
      ! set variables of mrm (redundant in case of coupling to mhm)
      call append(global_parameters, reshape(muskingumTravelTime_constant, (/1, nColPars/)))
      call append(global_parameters, reshape(muskingumTravelTime_riverLength, (/1, nColPars/)))
      call append(global_parameters, reshape(muskingumTravelTime_riverSlope, (/1, nColPars/)))
      call append(global_parameters, reshape(muskingumTravelTime_impervious, (/1, nColPars/)))
      call append(global_parameters, reshape(muskingumAttenuation_riverSlope, (/1, nColPars/)))

      call append(global_parameters_name, (/ &
              'muskingumTravelTime_constant   ', &
                      'muskingumTravelTime_riverLength', &
                      'muskingumTravelTime_riverSlope ', &
                      'muskingumTravelTime_impervious ', &
                      'muskingumAttenuation_riverSlope'/))
      ! adaptive timestep routing
    else if (processCase .eq. 2_i4) then
      processMatrix(8, 1) = processCase
      processMatrix(8, 2) = 1_i4
      processMatrix(8, 3) = processMatrix(8, 2)
      ! set variables of mrm (redundant in case of coupling to mhm)
      call append(global_parameters, reshape(streamflow_celerity, (/1, nColPars/)))
  

      call append(global_parameters_name, (/ &
              'streamflow_celerity'/))
     ! adaptive timestep routing - varying celerity
     else if (processCase .eq. 3_i4) then
        processMatrix(8, 1) = processCase
        processMatrix(8, 2) = 2_i4
        processMatrix(8, 3) = processMatrix(8, 2)
        ! set variables of mrm (redundant in case of coupling to mhm)
        call append(global_parameters, reshape(g1, (/1, nColPars/)))
        call append(global_parameters, reshape(g2, (/1, nColPars/)))

        call append(global_parameters_name, (/ &
             'g1', &
             'g2'/))
    end if
#endif

    ! check if parameter are in range
    if (.not. in_bound(global_parameters)) then
      call message('***ERROR: parameter in routing namelist out of bound in ', &
              trim(adjustl(file_namelist_param)))
      stop
    end if

    call close_nml(unamelist_param)

  end subroutine read_mrm_routing_params
end module mo_mrm_read_config
