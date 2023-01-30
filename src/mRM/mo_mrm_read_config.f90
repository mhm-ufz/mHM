!> \file mo_mrm_read_config.f90
!> \brief \copybrief mo_mrm_read_config
!> \details \copydetails mo_mrm_read_config

!> \brief read mRM config
!> \details This module contains all mRM subroutines related to reading the mRM configuration either from file or copy from mHM.
!> \authors Stephan Thober
!> \date Aug 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
module mo_mrm_read_config

  use mo_kind, only : i4, dp
  use mo_message, only: message, error_message

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

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober  Sep 2015 - removed stop condition when routing resolution is smaller than hydrologic resolution
  ! Stephan Thober  Oct 2015 - added NLoutputResults namelist, fileLatLon to directories_general namelist, and readLatLon flag
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_read_config(file_namelist, unamelist, file_namelist_param, unamelist_param, do_message)

    use mo_common_constants, only : maxNoDomains, nodata_i4
    use mo_common_mHM_mRM_read_config, only : common_check_resolution
    use mo_common_variables, only : ALMA_convention, domainMeta, processMatrix
    use mo_mrm_constants, only : maxNoGauges
    use mo_mrm_file, only : file_defOutput, udefOutput
    use mo_mrm_global_variables, only : InflowGauge, domainInfo_mRM, domain_mrm, &
                                        dirGauges, dirTotalRunoff, filenameTotalRunoff, dirBankfullRunoff, gauge, is_start, &
                                        nGaugesTotal, nGaugesLocal, nInflowGaugesTotal, outputFlxState_mrm, &
                                        timeStep_model_outputs_mrm, &
                                        varnameTotalRunoff, gw_coupling, &
                                        output_deflate_level_mrm, output_double_precision_mrm, output_time_reference_mrm, &
                                        readLatLon
    use mo_nml, only : close_nml, open_nml, position_nml
    use mo_string_utils, only : num2str

    implicit none

    character(*), intent(in) :: file_namelist, file_namelist_param

    integer, intent(in) :: unamelist, unamelist_param

    ! - flag for writing mHM standard messages
    logical, intent(in) :: do_message

    integer(i4), dimension(maxNoDomains) :: NoGauges_domain

    integer(i4), dimension(maxNoDomains, maxNoGauges) :: Gauge_id

    character(256), dimension(maxNoDomains, maxNoGauges) :: Gauge_filename

    integer(i4), dimension(maxNoDomains) :: NoInflowGauges_domain

    integer(i4), dimension(maxNoDomains, maxNoGauges) :: InflowGauge_id

    character(256), dimension(maxNoDomains, maxNoGauges) :: InflowGauge_filename

    logical, dimension(maxNoDomains, maxNoGauges) :: InflowGauge_Headwater

    integer(i4) :: domainID, iDomain

    integer(i4) :: iGauge

    integer(i4) :: idx

    character(256), dimension(maxNoDomains) :: dir_Gauges

    character(256), dimension(maxNoDomains) :: dir_Total_Runoff

    character(256), dimension(maxNoDomains) :: dir_Bankfull_Runoff

    logical :: file_exists

    type(domainInfo_mRM), pointer :: domain_mrm_iDomain


    ! namelist spatial & temporal resolution, optmization information
    namelist /mainconfig_mrm/ ALMA_convention, &
      filenameTotalRunoff, varnameTotalRunoff, gw_coupling
    ! namelist directories
    namelist /directories_mRM/ dir_Gauges, dir_Total_Runoff, dir_Bankfull_Runoff
    namelist /evaluation_gauges/ nGaugesTotal, NoGauges_domain, Gauge_id, gauge_filename
    ! namelist for inflow gauges
    namelist /inflow_gauges/ nInflowGaugesTotal, NoInflowGauges_domain, InflowGauge_id, &
            InflowGauge_filename, InflowGauge_Headwater
    ! name list regarding output
    namelist /NLoutputResults/ &
            output_deflate_level_mrm, &
            output_double_precision_mrm, &
            output_time_reference_mrm, &
            timeStep_model_outputs_mrm, &
            outputFlxState_mrm

    !===============================================================
    ! INITIALIZATION
    !===============================================================
    is_start = .True.
    nGaugesTotal = nodata_i4
    NoGauges_domain = nodata_i4
    Gauge_id = nodata_i4
    gauge_filename = num2str(nodata_i4)

    ! default arguments
    ALMA_convention = .false.
    filenameTotalRunoff = 'total_runoff'
    varnameTotalRunoff = 'total_runoff'
    gw_coupling = .false.

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

    allocate(dirGauges(domainMeta%nDomains), dirTotalRunoff(domainMeta%nDomains), dirBankfullRunoff(domainMeta%nDomains))
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      dirGauges(iDomain)         = dir_Gauges(domainID)
      dirTotalRunoff(iDomain)    = dir_Total_Runoff(domainID)
      dirBankfullRunoff(iDomain) = dir_Bankfull_Runoff(domainID)
    end do

    !===============================================================
    ! READ EVALUATION GAUGES
    !===============================================================
    call position_nml('evaluation_gauges', unamelist)
    read(unamelist, nml = evaluation_gauges)

    if (nGaugesTotal .GT. maxNoGauges) then
      call error_message('***ERROR: ', trim(file_namelist), ': Total number of evaluation gauges is restricted to', &
              num2str(maxNoGauges), raise=.false.)
      call error_message('          Error occured in namlist: evaluation_gauges')
    end if

    ! ToDo: check
    nGaugesLocal = 0
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      nGaugesLocal = nGaugesLocal + NoGauges_domain(domainID)
    end do
    ! End ToDo

    allocate(gauge%gaugeId(nGaugesLocal)) ; gauge%gaugeId = nodata_i4
    allocate(gauge%domainId(nGaugesLocal)) ; gauge%domainId = nodata_i4
    allocate(gauge%fName  (nGaugesLocal))
    if (nGaugesLocal > 0) then
      gauge%fName(1) = num2str(nodata_i4)
    end if
    allocate(domain_mrm(domainMeta%nDomains))

    idx = 0
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      domain_mrm_iDomain => domain_mrm(iDomain)
      ! initialize
      domain_mrm_iDomain%nGauges = nodata_i4
      allocate(domain_mrm_iDomain%gaugeIdList(maxval(NoGauges_domain(:))))
      domain_mrm_iDomain%gaugeIdList = nodata_i4
      allocate(domain_mrm_iDomain%gaugeIndexList(maxval(NoGauges_domain(:))))
      domain_mrm_iDomain%gaugeIndexList = nodata_i4
      allocate(domain_mrm_iDomain%gaugeNodeList(maxval(NoGauges_domain(:))))
      domain_mrm_iDomain%gaugeNodeList = nodata_i4
      ! check if NoGauges_domain has a valid value
      if (NoGauges_domain(domainID) .EQ. nodata_i4) then
        call error_message('***ERROR: ', trim(file_namelist), ': Number of evaluation gauges for subdomain ', &
                trim(adjustl(num2str(domainID))), ' is not defined!', raise=.false.)
        call error_message('          Error occured in namelist: evaluation_gauges')
      end if

      domain_mrm_iDomain%nGauges = NoGauges_domain(domainID)

      do iGauge = 1, NoGauges_domain(domainID)
        ! check if NoGauges_domain has a valid value
        if (Gauge_id(domainID, iGauge) .EQ. nodata_i4) then
          call error_message('***ERROR: ', trim(file_namelist), ': ID ', &
                  trim(adjustl(num2str(Gauge_id(domainID, iGauge)))), ' of evaluation gauge ', &
                  trim(adjustl(num2str(iGauge))), ' for subdomain ', &
                  trim(adjustl(num2str(iDomain))), ' is not defined!', raise=.false.)
          call error_message('          Error occured in namelist: evaluation_gauges')
        else if (trim(gauge_filename(domainID, iGauge)) .EQ. trim(num2str(nodata_i4))) then
          call error_message('***ERROR: ', trim(file_namelist), ': Filename of evaluation gauge ', &
                  trim(adjustl(num2str(iGauge))), ' for subdomain ', &
                  trim(adjustl(num2str(iDomain))), ' is not defined!', raise=.false.)
          call error_message('          Error occured in namelist: evaluation_gauges')
        end if
        !
        idx = idx + 1
        gauge%domainId(idx) = iDomain
        gauge%gaugeId(idx) = Gauge_id(domainID, iGauge)
        gauge%fname(idx) = trim(dirGauges(iDomain)) // trim(gauge_filename(domainID, iGauge))
        domain_mrm_iDomain%gaugeIdList(iGauge) = Gauge_id(domainID, iGauge)
        domain_mrm_iDomain%gaugeIndexList(iGauge) = idx
      end do
    end do

    if (nGaugesLocal .NE. idx) then
      call error_message('***ERROR: ', trim(file_namelist), ': Total number of evaluation gauges (', &
              trim(adjustl(num2str(nGaugesLocal))), &
              ') different from sum of gauges in subdomains (', trim(adjustl(num2str(idx))), ')!', raise=.false.)
      call error_message('          Error occured in namelist: evaluation_gauges')
    end if

    !===============================================================
    ! Read inflow gauge information
    !===============================================================

    nInflowGaugesTotal = 0
    NoInflowGauges_domain = 0
    InflowGauge_id = nodata_i4
    InflowGauge_filename = num2str(nodata_i4)

    call position_nml('inflow_gauges', unamelist)
    read(unamelist, nml = inflow_gauges)

    if (nInflowGaugesTotal .GT. maxNoGauges) then
      call error_message('***ERROR: ', trim(file_namelist), &
              ':read_gauge_lut: Total number of inflow gauges is restricted to', num2str(maxNoGauges), raise=.false.)
      call error_message('          Error occured in namlist: inflow_gauges')
    end if

    ! allocation - max() to avoid allocation with zero, needed for mhm call
    allocate(InflowGauge%gaugeId (max(1, nInflowGaugesTotal)))
    allocate(InflowGauge%domainId (max(1, nInflowGaugesTotal)))
    allocate(InflowGauge%fName   (max(1, nInflowGaugesTotal)))
    ! dummy initialization
    InflowGauge%gaugeId = nodata_i4
    InflowGauge%domainId = nodata_i4
    InflowGauge%fName = num2str(nodata_i4)

    idx = 0
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      domain_mrm_iDomain => domain_mrm(iDomain)

      allocate(domain_mrm_iDomain%InflowGaugeIdList    (max(1, maxval(NoInflowGauges_domain(:)))))
      allocate(domain_mrm_iDomain%InflowGaugeHeadwater (max(1, maxval(NoInflowGauges_domain(:)))))
      allocate(domain_mrm_iDomain%InflowGaugeIndexList (max(1, maxval(NoInflowGauges_domain(:)))))
      allocate(domain_mrm_iDomain%InflowGaugeNodeList  (max(1, maxval(NoInflowGauges_domain(:)))))
      ! dummy initialization
      domain_mrm_iDomain%nInflowGauges = 0
      domain_mrm_iDomain%InflowGaugeIdList = nodata_i4
      domain_mrm_iDomain%InflowGaugeHeadwater = .FALSE.
      domain_mrm_iDomain%InflowGaugeIndexList = nodata_i4
      domain_mrm_iDomain%InflowGaugeNodeList = nodata_i4
      ! no inflow gauge for subdomain i
      if (NoInflowGauges_domain(domainID) .EQ. nodata_i4) then
        NoInflowGauges_domain(domainID) = 0
      end if

      domain_mrm_iDomain%nInflowGauges = NoInflowGauges_domain(domainID)

      do iGauge = 1, NoInflowGauges_domain(domainID)
        ! check if NoInflowGauges_domain has a valid value
        if (InflowGauge_id(domainID, iGauge) .EQ. nodata_i4) then
          call error_message('***ERROR: ', trim(file_namelist), ':ID of inflow gauge ', &
                  trim(adjustl(num2str(iGauge))), ' for subdomain ', &
                  trim(adjustl(num2str(iDomain))), ' is not defined!', raise=.false.)
          call error_message('          Error occured in namlist: inflow_gauges')
        else if (trim(InflowGauge_filename(domainID, iGauge)) .EQ. trim(num2str(nodata_i4))) then
          call error_message('***ERROR: ', trim(file_namelist), ':Filename of inflow gauge ', &
                  trim(adjustl(num2str(iGauge))), ' for subdomain ', &
                  trim(adjustl(num2str(iDomain))), ' is not defined!', raise=.false.)
          call error_message('          Error occured in namlist: inflow_gauges')
        end if
        !
        idx = idx + 1
        InflowGauge%domainId(idx) = iDomain
        InflowGauge%gaugeId(idx) = InflowGauge_id(domainID, iGauge)
        InflowGauge%fname(idx) = trim(dirGauges(domainID)) // trim(InflowGauge_filename(domainID, iGauge))
        domain_mrm_iDomain%InflowGaugeIdList(iGauge) = InflowGauge_id(domainID, iGauge)
        domain_mrm_iDomain%InflowGaugeHeadwater(iGauge) = InflowGauge_Headwater(domainID, iGauge)
        domain_mrm_iDomain%InflowGaugeIndexList(iGauge) = idx
      end do
    end do

    if (nInflowGaugesTotal .NE. idx) then
      call error_message('***ERROR: ', trim(file_namelist), ': Total number of inflow gauges (', &
              trim(adjustl(num2str(nInflowGaugesTotal))), &
              ') different from sum of inflow gauges in subdomains (', trim(adjustl(num2str(idx))), ')!', raise=.false.)
      call error_message('          Error occured in namlist: inflow_gauges')
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
    output_deflate_level_mrm = 6
    output_time_reference_mrm = 0
    output_double_precision_mrm = .true.
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
      call message('    NetCDF deflate level: ', adjustl(trim(num2str(output_deflate_level_mrm))))
      if ( output_double_precision_mrm ) then
        call message('    NetCDF output precision: double')
      else
        call message('    NetCDF output precision: single')
      end if
      select case(output_time_reference_mrm)
        case(0)
          call message('    NetCDF output time reference point: start of time interval')
        case(1)
          call message('    NetCDF output time reference point: center of time interval')
        case(2)
          call message('    NetCDF output time reference point: end of time interval')
      end select
      call message('    FLUXES:')
      if (outputFlxState_mrm(1)) then
        call message('      routed streamflow      (L11_qMod)                [m3 s-1]')
      end if
      if (outputFlxState_mrm(2)) then
        call message('      river temperature      (RivTemp)                 [deg C]')
      end if
      if (gw_coupling) then
        call message('      river head             (river_head)              [m]')
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
    use mo_nml, only : close_nml, open_nml, position_nml

    implicit none

    ! it is the default case, should be one
    integer(i4), intent(in) :: processCase

    ! file name containing parameter namelist
    character(*), intent(in) :: file_namelist_param

    ! file name id containing parameter namelist
    integer(i4), intent(in) :: unamelist_param

    ! equals sum of previous parameters
    integer(i4) :: start_index

    real(dp), dimension(nColPars) :: muskingumTravelTime_constant

    real(dp), dimension(nColPars) :: muskingumTravelTime_riverLength

    real(dp), dimension(nColPars) :: muskingumTravelTime_riverSlope

    real(dp), dimension(nColPars) :: muskingumTravelTime_impervious

    real(dp), dimension(nColPars) :: muskingumAttenuation_riverSlope

    real(dp), dimension(nColPars) :: streamflow_celerity
    real(dp), dimension(nColPars) :: slope_factor

    namelist /routing1/ muskingumTravelTime_constant, muskingumTravelTime_riverLength, &
            muskingumTravelTime_riverSlope, muskingumTravelTime_impervious, muskingumAttenuation_riverSlope
    namelist /routing2/ streamflow_celerity
    namelist /routing3/ slope_factor
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
        processMatrix(8, 2) = 1_i4
        processMatrix(8, 3) = sum(processMatrix(1:8, 2))
        start_index         = processMatrix(8, 3) - processMatrix(8, 2)
        global_parameters(start_index + 1, :) = slope_factor

        global_parameters_name(start_index + 1 : start_index + processMatrix(8,2)) = (/ &
             'slope_factor'/)
    end if

    ! check if parameter are in range
    if (.not. in_bound(global_parameters)) then
      call error_message('***ERROR: parameter in routing namelist out of bound in ', &
              trim(adjustl(file_namelist_param)))
    end if

    call close_nml(unamelist_param)

  end subroutine read_mrm_routing_params
end module mo_mrm_read_config
