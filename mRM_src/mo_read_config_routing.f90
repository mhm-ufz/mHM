module mo_read_config_routing
  use mo_kind, only : i4, dp
  implicit none
  public :: read_config_routing
  public :: read_routing_params
contains
  subroutine read_config_routing(do_routing)
    use mo_message,          only: message
    use mo_nml,              only: position_nml
    use mo_mrm_constants, only : nodata_i4, &
         maxNoGauges    ! maximum number of allowed gauges
    use mo_file,          only : unamelist
    use mo_string_utils,  only : num2str
    use mo_global_variables_routing, only : nGaugesTotal, gauge, & ! number of evaluation gauges and gauge informations 
         nInflowGaugesTotal, InflowGauge, & ! number of inflow gauges and gauge informations
         dirGauges, & ! directory of gauge files
         is_start ! flag for first timestep
    !ST The following dependency has to be changed
    use mo_global_variables, only: basin, nbasins
    use mo_mhm_constants, only: maxNoBasins ! maximum number of allowed basins

    
    implicit none
    ! input variables
    logical, intent(in) :: do_routing
    ! local variables
    integer(i4),    dimension(maxNoBasins)             :: NoGauges_basin
    integer(i4),    dimension(maxNoBasins,maxNoGauges) :: Gauge_id
    character(256), dimension(maxNoGauges,maxNoGauges) :: Gauge_filename
    integer(i4),    dimension(maxNoBasins)             :: NoInflowGauges_basin
    integer(i4),    dimension(maxNoBasins,maxNoGauges) :: InflowGauge_id
    character(256), dimension(maxNoGauges,maxNoGauges) :: InflowGauge_filename
    logical,        dimension(maxNoBasins,maxNoGauges) :: InflowGauge_Headwater
    integer(i4)                                        :: i_basin
    integer(i4)                                        :: i_gauge
    integer(i4)                                        :: idx
    
    !
    ! namelist for evaluation gauges
    namelist /evaluation_gauges/ nGaugesTotal, NoGauges_basin, Gauge_id, gauge_filename
    ! namelist for inflow gauges
    namelist /inflow_gauges/ nInflowGaugesTotal, NoInflowGauges_basin, InflowGauge_id, &
         InflowGauge_filename, InflowGauge_Headwater
    ! initialization
    if (do_routing) then
       is_start = .True.
       nGaugesTotal   = nodata_i4
       NoGauges_basin = nodata_i4
       Gauge_id       = nodata_i4
       gauge_filename = num2str(nodata_i4)

       call position_nml('evaluation_gauges', unamelist)
       read(unamelist, nml=evaluation_gauges)

       if (nGaugesTotal .GT. maxNoGauges) then
          call message()
          call message('***ERROR: mhm.nml: Total number of evaluation gauges is restricted to', num2str(maxNoGauges))
          call message('          Error occured in namlist: evaluation_gauges')
          stop
       end if

       allocate(gauge%gaugeId        (nGaugesTotal))                       ; gauge%gaugeId        = nodata_i4
       allocate(gauge%basinId        (nGaugesTotal))                       ; gauge%basinId        = nodata_i4
       allocate(gauge%fName          (nGaugesTotal))                       ; gauge%fName(1)       = num2str(nodata_i4)
       allocate(basin%nGauges        (nBasins                           )) ; basin%nGauges        = nodata_i4
       allocate(basin%gaugeIdList    (nBasins, maxval(NoGauges_basin(:)))) ; basin%gaugeIdList    = nodata_i4
       allocate(basin%gaugeIndexList (nBasins, maxval(NoGauges_basin(:)))) ; basin%gaugeIndexList = nodata_i4
       allocate(basin%gaugeNodeList  (nBasins, maxval(NoGauges_basin(:)))) ; basin%gaugeNodeList  = nodata_i4

       idx = 0
       do i_basin = 1, nBasins
          ! check if NoGauges_basin has a valid value
          if ( NoGauges_basin(i_basin) .EQ. nodata_i4 ) then
             call message()
             call message('***ERROR: mhm.nml: Number of evaluation gauges for subbasin ', &
                  trim(adjustl(num2str(i_basin))),' is not defined!')
             call message('          Error occured in namlist: evaluation_gauges')
             stop
          end if

          basin%nGauges(i_basin)          = NoGauges_basin(i_basin)

          do i_gauge = 1, NoGauges_basin(i_basin)
             ! check if NoGauges_basin has a valid value
             if (Gauge_id(i_basin,i_gauge) .EQ. nodata_i4) then
                call message()
                call message('***ERROR: mhm.nml: ID of evaluation gauge ',        &
                     trim(adjustl(num2str(i_gauge))),' for subbasin ', &
                     trim(adjustl(num2str(i_basin))),' is not defined!')
                call message('          Error occured in namlist: evaluation_gauges')
                stop
             else if (trim(gauge_filename(i_basin,i_gauge)) .EQ. trim(num2str(nodata_i4))) then
                call message()
                call message('***ERROR: mhm.nml: Filename of evaluation gauge ', &
                     trim(adjustl(num2str(i_gauge))),' for subbasin ',  &
                     trim(adjustl(num2str(i_basin))),' is not defined!')
                call message('          Error occured in namlist: evaluation_gauges')
                stop
             end if
             !
             idx = idx + 1
             gauge%basinId(idx)                    = i_basin
             gauge%gaugeId(idx)                    = Gauge_id(i_basin,i_gauge)
             gauge%fname(idx)                      = trim(dirGauges(i_basin)) // trim(gauge_filename(i_basin,i_gauge))
             basin%gaugeIdList(i_basin,i_gauge)    = Gauge_id(i_basin,i_gauge)
             basin%gaugeIndexList(i_basin,i_gauge) = idx
          end do
       end do

       if ( nGaugesTotal .NE. idx) then
          call message()
          call message('***ERROR: mhm.nml: Total number of evaluation gauges (', trim(adjustl(num2str(nGaugesTotal))), &
               ') different from sum of gauges in subbasins (', trim(adjustl(num2str(idx))), ')!')
          call message('          Error occured in namlist: evaluation_gauges')
          stop
       end if

    end if

    !===============================================================
    ! Read inflow gauge information
    !===============================================================

    nInflowGaugesTotal   = 0
    NoInflowGauges_basin = 0
    InflowGauge_id       = nodata_i4
    InflowGauge_filename = num2str(nodata_i4)

    if (do_routing) then
       call position_nml('inflow_gauges', unamelist)
       read(unamelist, nml=inflow_gauges)

       if (nInflowGaugesTotal .GT. maxNoGauges) then
          call message()
          call message('***ERROR: mhm.nml:read_gauge_lut: Total number of inflow gauges is restricted to', num2str(maxNoGauges))
          call message('          Error occured in namlist: inflow_gauges')
          stop
       end if
    end if

    ! allocation - max() to avoid allocation with zero, needed for mhm call
    allocate(InflowGauge%gaugeId        (max(1,nInflowGaugesTotal)))
    allocate(InflowGauge%basinId        (max(1,nInflowGaugesTotal)))
    allocate(InflowGauge%fName          (max(1,nInflowGaugesTotal)))
    allocate(basin%nInflowGauges        (nBasins                                 ))
    allocate(basin%InflowGaugeIdList    (nBasins, max(1, maxval(NoInflowGauges_basin(:)))))
    allocate(basin%InflowGaugeHeadwater (nBasins, max(1, maxval(NoInflowGauges_basin(:)))))
    allocate(basin%InflowGaugeIndexList (nBasins, max(1, maxval(NoInflowGauges_basin(:)))))
    allocate(basin%InflowGaugeNodeList  (nBasins, max(1, maxval(NoInflowGauges_basin(:)))))
    ! dummy initialization
    InflowGauge%gaugeId        = nodata_i4
    InflowGauge%basinId        = nodata_i4
    InflowGauge%fName          = num2str(nodata_i4)
    basin%nInflowGauges        = 0
    basin%InflowGaugeIdList    = nodata_i4
    basin%InflowGaugeHeadwater = .FALSE.
    basin%InflowGaugeIndexList = nodata_i4
    basin%InflowGaugeNodeList  = nodata_i4

    if (do_routing) then
       idx = 0
       do i_basin = 1, nBasins

          ! no inflow gauge for subbasin i
          if (NoInflowGauges_basin(i_basin) .EQ. nodata_i4) then
             NoInflowGauges_basin(i_basin)       = 0
          end if

          basin%nInflowGauges(i_basin) = NoInflowGauges_basin(i_basin)

          do i_gauge = 1, NoInflowGauges_basin(i_basin)
             ! check if NoInflowGauges_basin has a valid value
             if (InflowGauge_id(i_basin,i_gauge) .EQ. nodata_i4) then
                call message()
                call message('***ERROR: mhm.nml:ID of inflow gauge ',        &
                     trim(adjustl(num2str(i_gauge))),' for subbasin ', &
                     trim(adjustl(num2str(i_basin))),' is not defined!')
                call message('          Error occured in namlist: inflow_gauges')
                stop
             else if (trim(InflowGauge_filename(i_basin,i_gauge)) .EQ. trim(num2str(nodata_i4))) then
                call message()
                call message('***ERROR: mhm.nml:Filename of inflow gauge ', &
                     trim(adjustl(num2str(i_gauge))),' for subbasin ',  &
                     trim(adjustl(num2str(i_basin))),' is not defined!')
                call message('          Error occured in namlist: inflow_gauges')
                stop
             end if
             !
             idx = idx + 1
             InflowGauge%basinId(idx)                    = i_basin
             InflowGauge%gaugeId(idx)                    = InflowGauge_id(i_basin,i_gauge)
             InflowGauge%fname(idx)                      = trim(dirGauges(i_basin)) // trim(InflowGauge_filename(i_basin,i_gauge))
             basin%InflowGaugeIdList(i_basin,i_gauge)    = InflowGauge_id(i_basin,i_gauge)
             basin%InflowGaugeHeadwater(i_basin,i_gauge) = InflowGauge_Headwater(i_basin,i_gauge)
             basin%InflowGaugeIndexList(i_basin,i_gauge) = idx
          end do
       end do

       if ( nInflowGaugesTotal .NE. idx) then
          call message()
          call message('***ERROR: mhm.nml: Total number of inflow gauges (', trim(adjustl(num2str(nInflowGaugesTotal))), &
               ') different from sum of inflow gauges in subbasins (', trim(adjustl(num2str(idx))), ')!')
          call message('          Error occured in namlist: inflow_gauges')
          stop
       end if
    end if

  end subroutine read_config_routing

  subroutine read_routing_params(processCase)
    use mo_file, only: file_namelist_param, unamelist_param ! file containing parameter values
    use mo_append, only: append
    use mo_nml, only: position_nml
    use mo_message, only: message
    use mo_mhm_constants, only: nColPars ! number of properties of the global variables
    !ST the following dependency has to be removed
    use mo_global_variables, only :  processMatrix, & ! process configuration
         global_parameters, & ! global parameters
         global_parameters_name ! clear names of global parameters
    implicit none
    ! input variables
    integer(i4), intent(in) :: processCase ! it is the default case should be one
    ! local variables
    real(dp), dimension(nColPars) :: muskingumTravelTime_constant
    real(dp), dimension(nColPars) :: muskingumTravelTime_riverLength
    real(dp), dimension(nColPars) :: muskingumTravelTime_riverSlope
    real(dp), dimension(nColPars) :: muskingumTravelTime_impervious
    real(dp), dimension(nColPars) :: muskingumAttenuation_riverSlope

    namelist /routing1/ muskingumTravelTime_constant, muskingumTravelTime_riverLength, &
         muskingumTravelTime_riverSlope, muskingumTravelTime_impervious, muskingumAttenuation_riverSlope
    !
    call position_nml('routing1', unamelist_param)
    read(unamelist_param, nml=routing1)
    processMatrix(8, 1) = processCase
    processMatrix(8, 2) = 5_i4
    processMatrix(8, 3) = sum(processMatrix(1:8, 2))
    call append(global_parameters, reshape(muskingumTravelTime_constant,    (/1, nColPars/)))
    call append(global_parameters, reshape(muskingumTravelTime_riverLength, (/1, nColPars/)))
    call append(global_parameters, reshape(muskingumTravelTime_riverSlope,  (/1, nColPars/)))
    call append(global_parameters, reshape(muskingumTravelTime_impervious,  (/1, nColPars/)))
    call append(global_parameters, reshape(muskingumAttenuation_riverSlope, (/1, nColPars/)))

    call append(global_parameters_name, (/ &
         'muskingumTravelTime_constant   ', &
         'muskingumTravelTime_riverLength', &
         'muskingumTravelTime_riverSlope ', &
         'muskingumTravelTime_impervious ', &
         'muskingumAttenuation_riverSlope'/))

    ! check if parameter are in range
    if ( .not. in_bound(global_parameters) ) then
       call message('***ERROR: parameter in namelist "routing1" out of bound in ', &
            trim(adjustl(file_namelist_param)))
       stop
    end if

  end subroutine read_routing_params

  ! --------------------------------------------------------------------------------
  ! private funtions and subroutines, DUPLICATED FROM mo_read_config.f90
  ! --------------------------------------------------------------------------------

  function in_bound(params)
    real(dp), dimension(:,:), intent(in) :: params ! parameter:
    !                                              !   col_1=Lower bound,
    !                                              !   col_2=Upper bound
    !                                              !   col_3=initial
    logical :: in_bound

    if ( any(params(:,3) .lt. params(:,1)) .or. any(params(:,3) .gt. params(:,2)) ) then
       in_bound=.false.
    else
       in_bound=.true.
    end if

  end function in_bound

end module mo_read_config_routing
