!>       \file mo_mrm_init.f90

!>       \brief Wrapper for initializing Routing.

!>       \details Calling all routines to initialize all mRM variables

!>       \authors Luis Samaniego, Rohini Kumar and Stephan Thober

!>       \date Aug 2015

! Modifications:

MODULE mo_mrm_init

    use mo_common_variables, only : dirOut

  ! This module sets the river network characteristics and routing order.

  ! Written  Luis Samaniego, Mar 2005

  IMPLICIT NONE

  public :: mrm_init
  public :: variables_default_init_routing

  private

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_init

  !    PURPOSE
  !>       \brief Initialize all mRM variables at all levels (i.e., L0, L1, and L11).

  !>       \details Initialize all mRM variables at all levels (i.e., L0, L1, and L11)
  !>       either with default values or with values from restart file. The L0 mask (L0_mask),
  !>       L0 elevation (L0_elev), and L0 land cover (L0_LCover) can be provided as optional
  !>       variables to save memory because these variable will then not be read in again.

  !    INTENT(IN)
  !>       \param[in] "character(*) :: file_namelist, file_namelist_param"
  !>       \param[in] "integer :: unamelist, unamelist_param"
  !>       \param[in] "character(*) :: file_namelist, file_namelist_param"
  !>       \param[in] "integer :: unamelist, unamelist_param"

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Stephan Thober Sep 2015 - added L0_mask, L0_elev, and L0_LCover
  ! Stephan Thober May 2016 - added warning message in case no gauge is found in modelling domain
  ! Matthias Kelbling Aug 2017 - added L11_flow_accumulation to Initialize Stream Netwo  
  ! Lennart Schueler May 2018 - added initialization for groundwater coupling
  ! Stephan Thober Jun 2018 - refactored for mpr_extract version


  subroutine mrm_init(file_namelist, unamelist, file_namelist_param, unamelist_param)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_mHM_mRM_read_config, only : check_optimization_settings, common_mHM_mRM_read_config
    use mo_common_mHM_mRM_variables, only : dirRestartIn, mrm_coupling_mode, read_restart, &
                                            resolutionRouting
    use mo_common_read_config, only : common_read_config
    use mo_common_restart, only : read_grid_info
    use mo_common_variables, only : L0_Basin, global_parameters, l0_l1_remap, level0, level1, nBasins, &
                                    processMatrix, resolutionHydrology
    use mo_grid, only : L0_grid_setup, init_lowres_level, set_basin_indices
    use mo_kind, only : i4
    use mo_message, only : message
    use mo_mrm_global_variables, only : basin_mrm, &
                                        l0_l11_remap, l1_l11_remap, level11, gw_coupling, &
                                        L0_river_head_mon_sum
    use mo_mrm_net_startup, only : L11_flow_direction, L11_fraction_sealed_floodplain, &
                                   L11_link_location, L11_routing_order, L11_set_drain_outlet_gauges, &
                                   L11_set_network_topology, L11_stream_features, l11_l1_mapping
    use mo_mrm_read_config, only : mrm_read_config
    use mo_mrm_read_data, only : mrm_read_L0_data, mrm_read_discharge, &
                                 mrm_read_total_runoff, mrm_read_bankfull_runoff
    use mo_mrm_restart, only : mrm_read_restart_config
    use mo_read_latlon, only : read_latlon
    use mo_mrm_river_head, only: init_masked_zeros_l0, create_output, calc_channel_elevation
    use mo_mrm_mpr, only : mrm_init_param

    implicit none

    character(*), intent(in) :: file_namelist, file_namelist_param

    integer, intent(in) :: unamelist, unamelist_param

    ! start and end index for routing parameters
    integer(i4) :: iStart, iEnd

    integer(i4) :: iBasin, gauge_counter

    logical :: ReadLatLon


    if (mrm_coupling_mode .eq. 0_i4) then
      call common_read_config(file_namelist, unamelist)
      call common_mHM_mRM_read_config(file_namelist, unamelist)
      !-----------------------------------------------------------
      ! PRINT STARTUP MESSAGE
      !-----------------------------------------------------------
      call print_startup_message(file_namelist, file_namelist_param)
    else
      call message('')
      call message('  Inititalize mRM')
    end if

    ! read config for mrm, readlatlon is set here depending on whether output is needed
    call mrm_read_config(file_namelist, unamelist, file_namelist_param, unamelist_param, &
            (mrm_coupling_mode .eq. 0_i4), ReadLatLon)

    ! this was moved here, because it depends on global_parameters that are only set in mrm_read_config
    if (mrm_coupling_mode .eq. 0_i4) then
      call check_optimization_settings()
      allocate(l0_l1_remap(nBasins))
      allocate(level1(nBasins))
      !-----------------------------------------------------------
      ! CONFIG OUTPUT
      !-----------------------------------------------------------
      call config_output()
    end if

    ! ----------------------------------------------------------
    ! READ DATA
    ! ----------------------------------------------------------
    allocate(level11(nBasins))
    allocate(l0_l11_remap(nBasins))
    allocate(l1_l11_remap(nBasins))

    if (.not. read_restart) then
      ! read all (still) necessary level 0 data
      if (processMatrix(8, 1) .eq. 1_i4) call mrm_read_L0_data(mrm_coupling_mode .eq. 0_i4, ReadLatLon, .true.)
      if (processMatrix(8, 1) .eq. 2_i4) call mrm_read_L0_data(mrm_coupling_mode .eq. 0_i4, ReadLatLon, .false.)
      if (processMatrix(8, 1) .eq. 3_i4) call mrm_read_L0_data(mrm_coupling_mode .eq. 0_i4, ReadLatLon, .false.)
    end if

    do iBasin = 1, nBasins
      if (read_restart) then
        ! this reads the basin properties
        if (mrm_coupling_mode .eq. 0_i4) then
          call read_grid_info(iBasin, dirRestartIn(iBasin), "1", "mRM", level1(iBasin))
        end if
        call read_grid_info(iBasin, dirRestartIn(iBasin), "11", "mRM", level11(iBasin))
        call mrm_read_restart_config(iBasin, dirRestartIn(iBasin))
      else
        if (iBasin .eq. 1) then
          call L0_check_input_routing(L0_Basin(iBasin))
          if (mrm_coupling_mode .eq. 0_i4) then
            call L0_grid_setup(level0(L0_Basin(iBasin)))
          end if
        else if ((L0_Basin(iBasin) .ne. L0_Basin(iBasin - 1))) then
          call L0_check_input_routing(L0_Basin(iBasin))
          if (mrm_coupling_mode .eq. 0_i4) then
            call L0_grid_setup(level0(L0_Basin(iBasin)))
          end if
        end if

        if (mrm_coupling_mode .eq. 0_i4) then
          call init_lowres_level(level0(L0_Basin(iBasin)), resolutionHydrology(iBasin), &
                  level1(iBasin), l0_l1_remap(iBasin))
        end if
        call init_lowres_level(level0(L0_Basin(iBasin)), resolutionRouting(iBasin), &
                level11(iBasin), l0_l11_remap(iBasin))
        call init_lowres_level(level1(iBasin), resolutionRouting(iBasin), &
                level11(iBasin), l1_l11_remap(iBasin))
        call L11_L1_mapping(iBasin)

        if (ReadLatLon) then
          ! read lat lon coordinates of each basin
          call read_latlon(iBasin, "lon_l11", "lat_l11", "level11", level11(iBasin))
        else
          ! allocate the memory and set to nodata
          allocate(level11(iBasin)%x(level11(iBasin)%nrows, level11(iBasin)%ncols))
          allocate(level11(iBasin)%y(level11(iBasin)%nrows, level11(iBasin)%ncols))
          level11(iBasin)%x = nodata_dp
          level11(iBasin)%y = nodata_dp
        end if
      end if
    end do

    call set_basin_indices(level11)
    call set_basin_indices(level1)

    ! ----------------------------------------------------------
    ! INITIALIZE STATES AND AUXILLIARY VARIABLES
    ! ----------------------------------------------------------
    do iBasin = 1, nBasins
      call variables_alloc_routing(iBasin)
    end do

    ! ----------------------------------------------------------
    ! INITIALIZE STREAM NETWORK
    ! ----------------------------------------------------------
    do iBasin = 1, nBasins
      if (.not. read_restart) then
        call L11_flow_direction(iBasin)
        call L11_set_network_topology(iBasin)
        call L11_routing_order(iBasin)
        call L11_link_location(iBasin)
        call L11_set_drain_outlet_gauges(iBasin)
        ! stream characteristics
        call L11_stream_features(iBasin)
      end if
    end do

    ! ----------------------------------------------------------
    ! INITIALIZE PARAMETERS
    ! ----------------------------------------------------------
    do iBasin = 1, nBasins
      iStart = processMatrix(8, 3) - processMatrix(8, 2) + 1
      iEnd = processMatrix(8, 3)
      call mrm_init_param(iBasin, global_parameters(iStart : iEnd, 3))
    end do

    ! check whether there are gauges within the modelling domain
    if (allocated(basin_mrm)) then
      gauge_counter = 0
      do iBasin = 1, nBasins
        if (.not. all(basin_mrm(iBasin)%gaugeNodeList .eq. nodata_i4)) then
          gauge_counter = gauge_counter + 1
        end if
      end do
      if (gauge_counter .lt. 1) then
        call message('')
        call message('    WARNING: no gauge found within modelling domain')
      end if
    end if
    ! mpr-like definiton of sealed floodplain fraction
    if ((processMatrix(8, 1) .eq. 1_i4) .and. (.not. read_restart)) then
      call L11_fraction_sealed_floodplain(2_i4, .true.)
    else
      ! dummy initialization
      call L11_fraction_sealed_floodplain(2_i4, .false.)
    end if

    ! -------------------------------------------------------
    ! READ INPUT DATA AND OBSERVED DISCHARGE DATA
    ! -------------------------------------------------------
    ! read simulated runoff at level 1
    if (mrm_coupling_mode .eq. 0_i4) then
      do iBasin = 1, nBasins
        call mrm_read_total_runoff(iBasin)
      end do
    end if
    ! discharge data
    call mrm_read_discharge()

    if (gw_coupling) then
        do iBasin = 1, nBasins
            call init_masked_zeros_l0(iBasin, L0_river_head_mon_sum)
            call mrm_read_bankfull_runoff(iBasin)
            call create_output(iBasin, dirOut(iBasin))
        end do
        call calc_channel_elevation()
    end if

    call message('')
    call message('  Finished Initialization of mRM')

  end subroutine mrm_init

  
  !===============================================================
  ! PRINT STARTUP MESSAGE
  !===============================================================
  !    NAME
  !        print_startup_message

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "character(*) :: file_namelist, file_namelist_param"
  !>       \param[in] "character(*) :: file_namelist, file_namelist_param"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine print_startup_message(file_namelist, file_namelist_param)

    use mo_kind, only : i4
    use mo_message, only : message, message_text
    use mo_mrm_file, only : file_defOutput, file_main, version, version_date
    use mo_string_utils, only : num2str, separator

    implicit none

    character(*), intent(in) :: file_namelist, file_namelist_param

    ! Date and time
    integer(i4), dimension(8) :: datetime


    call message(separator)
    call message('              mRM-UFZ')
    call message()
    call message('    MULTISCALE ROUTING MODEL')
    call message('           Version ', trim(version))
    call message('           ', trim(version_date))
    call message()
    call message('Made available by S. Thober & M. Cuntz')
    call message()
    call message('Based on mHM-UFZ by L. Samaniego & R. Kumar')

    call message(separator)

    call message()
    call date_and_time(values = datetime)
    message_text = trim(num2str(datetime(3), '(I2.2)')) // "." // trim(num2str(datetime(2), '(I2.2)')) &
            // "." // trim(num2str(datetime(1), '(I4.4)')) // " " // trim(num2str(datetime(5), '(I2.2)')) &
            // ":" // trim(num2str(datetime(6), '(I2.2)')) // ":" // trim(num2str(datetime(7), '(I2.2)'))
    call message('Start at ', trim(message_text), '.')
    call message('Using main file ', trim(file_main), ' and namelists: ')
    call message('     ', trim(file_namelist))
    call message('     ', trim(file_namelist_param))
    call message('     ', trim(file_defOutput), ' (if it is given)')
    call message()

  end subroutine print_startup_message


  !---------------------------------------------------------------
  ! Config output
  !---------------------------------------------------------------
  !    NAME
  !        config_output

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine config_output

    use mo_common_variables, only : dirLCover, dirMorpho, dirOut, nBasins
    use mo_kind, only : i4
    use mo_message, only : message
    use mo_mrm_file, only : file_defOutput, file_namelist_mrm, file_namelist_param_mrm
    use mo_mrm_global_variables, only : basin_mrm, &
                                        dirGauges
    use mo_string_utils, only : num2str

    implicit none

    integer(i4) :: iBasin

    integer(i4) :: jj


    !
    call message()
    call message('Read namelist file: ', trim(file_namelist_mrm))
    call message('Read namelist file: ', trim(file_namelist_param_mrm))
    call message('Read namelist file: ', trim(file_defOutput), ' (if it is given)')

    call message()
    call message('  # of basins:         ', trim(num2str(nbasins)))
    call message()
    call message('  Input data directories:')
    do iBasin = 1, nbasins
      call message('  --------------')
      call message('      BASIN                   ', num2str(iBasin, '(I3)'))
      call message('  --------------')
      call message('    Morphological directory:    ', trim(dirMorpho(iBasin)))
      call message('    Land cover directory:       ', trim(dirLCover(iBasin)))
      call message('    Discharge directory:        ', trim(dirGauges(iBasin)))
      call message('    Output directory:           ', trim(dirOut(iBasin)))
      call message('    Evaluation gauge            ', 'ID')
      do jj = 1, basin_mrm(iBasin)%nGauges
        call message('    ', trim(adjustl(num2str(jj))), '                           ', &
                trim(adjustl(num2str(basin_mrm(iBasin)%gaugeIdList(jj)))))
      end do
      if (basin_mrm(iBasin)%nInflowGauges .GT. 0) then
        call message('    Inflow gauge              ', 'ID')
        do jj = 1, basin_mrm(iBasin)%nInflowGauges
          call message('    ', trim(adjustl(num2str(jj))), '                         ', &
                  trim(adjustl(num2str(basin_mrm(iBasin)%InflowGaugeIdList(jj)))))
        end do
      end if
    end do
  end subroutine config_output


  ! ------------------------------------------------------------------

  !    NAME
  !        variables_default_init_routing

  !    PURPOSE
  !>       \brief Default initalization mRM related L11 variables

  !>       \details Default initalization of mHM related L11 variables (e.g., states,
  !>       fluxes, and parameters) as per given constant values given in mo_mhm_constants.
  !>       Variables initalized here is defined in the mo_global_variables.f90 file.
  !>       Only Variables that are defined in the variables_alloc subroutine are
  !>       intialized here.
  !>       If a variable is added or removed here, then it also has to be added or removed
  !>       in the subroutine state_variables_set in the module mo_restart and in the
  !>       subroutine set_state in the module mo_set_netcdf_restart.
  !>       ADDITIONAL INFORMATION
  !>       variables_default_init_routing


  !>       call variables_default_init()
  !>       \authors  Stephan Thober, Rohini Kumar, and Juliane Mai
  !>       \date    Aug 2015

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine variables_default_init_routing

    use mo_common_constants, only : P1_InitStateFluxes
    use mo_mrm_global_variables, only : L11_C1, L11_C2, L11_K, L11_Qmod, L11_qOUT, L11_qTIN, L11_qTR, L11_xi

    implicit none


    !-------------------------------------------
    ! L11 ROUTING STATE VARIABLES, FLUXES AND
    !             PARAMETERS
    !-------------------------------------------
    ! simulated discharge at each node
    L11_Qmod = P1_InitStateFluxes
    ! Total outflow from cells L11 at time tt
    L11_qOUT = P1_InitStateFluxes
    ! Total discharge inputs at t-1 and t
    L11_qTIN = P1_InitStateFluxes
    !  Routed outflow leaving a node
    L11_qTR = P1_InitStateFluxes
    ! kappa: Muskingum travel time parameter.
    L11_K = P1_InitStateFluxes
    ! xi:    Muskingum diffusion parameter
    L11_xi = P1_InitStateFluxes
    ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
    L11_C1 = P1_InitStateFluxes
    ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
    L11_C2 = P1_InitStateFluxes
  end subroutine variables_default_init_routing


  ! --------------------------------------------------------------------------
  ! L0_check_input_routing
  ! --------------------------------------------------------------------------
  !    NAME
  !        L0_check_input_routing

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: L0Basin_iBasin"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine L0_check_input_routing(L0Basin_iBasin)

    use mo_common_constants, only : nodata_i4
    use mo_common_variables, only : level0
    use mo_kind, only : i4
    use mo_message, only : message, message_text
    use mo_mrm_global_variables, only : L0_fAcc, L0_fDir
    use mo_string_utils, only : num2str

    implicit none

    integer(i4), intent(in) :: L0Basin_iBasin

    integer(i4) :: k


    do k = level0(L0Basin_iBasin)%iStart, level0(L0Basin_iBasin)%iEnd
      ! flow direction [-]
      if (L0_fDir(k) .eq. nodata_i4) then
        message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(L0Basin_iBasin, '(I5)'))
        call message(' Error: flow direction has missing value within the valid masked area at cell in basin ', &
                trim(message_text))
        stop
      end if
      ! flow accumulation [-]
      if (L0_fAcc(k) .eq. nodata_i4) then
        message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(L0Basin_iBasin, '(I5)'))
        call message(' Error: flow accumulation has missing values within the valid masked area at cell in basin ', &
                trim(message_text))
        stop
      end if
    end do

  end subroutine L0_check_input_routing


  ! --------------------------------------------------------------------------
  ! L11 ROUTING STATE VARIABLES, FLUXES AND
  !             PARAMETERS
  ! --------------------------------------------------------------------------
  !    NAME
  !        variables_alloc_routing

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iBasin"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine variables_alloc_routing(iBasin)

    use mo_append, only : append
    use mo_kind, only : dp, i4
    use mo_mrm_constants, only : nRoutingStates
    use mo_common_variables, only : level0, L0_Basin
    use mo_mrm_global_variables, only : L11_C1, L11_C2, L11_K, &
         L11_Qmod, L11_qOUT, L11_qTIN, L11_qTR, L11_xi, &
         level11, L11_celerity, L0_celerity

    implicit none

    integer(i4), intent(in) :: iBasin

    real(dp), dimension(:), allocatable :: dummy_Vector11

    real(dp), dimension(:, :), allocatable :: dummy_Matrix11_IT


    ! dummy vector and matrix
    allocate(dummy_Vector11   (level11(iBasin)%nCells))
    allocate(dummy_Matrix11_IT(level11(iBasin)%nCells, nRoutingStates))

    ! simulated discharge at each node
    dummy_Vector11(:) = 0.0_dp
    call append(L11_Qmod, dummy_Vector11)

    ! Total outflow from cells L11 at time tt
    dummy_Vector11(:) = 0.0_dp
    call append(L11_qOUT, dummy_Vector11)

    ! Total discharge inputs at t-1 and t
    dummy_Matrix11_IT(:, :) = 0.0_dp
    call append(L11_qTIN, dummy_Matrix11_IT)

    !  Routed outflow leaving a node
    dummy_Matrix11_IT(:, :) = 0.0_dp
    call append(L11_qTR, dummy_Matrix11_IT)

    ! kappa: Muskingum travel time parameter.
    dummy_Vector11(:) = 0.0_dp
    call append(L11_K, dummy_Vector11)

    ! xi:    Muskingum diffusion parameter
    dummy_Vector11(:) = 0.0_dp
    call append(L11_xi, dummy_Vector11)

    ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
    dummy_Vector11(:) = 0.0_dp
    call append(L11_C1, dummy_Vector11)

    ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
    dummy_Vector11(:) = 0.0_dp
    call append(L11_C2, dummy_Vector11)

    ! Celerity at each link
    dummy_Vector11(:) = 0.0_dp
    call append(L11_celerity, dummy_Vector11)

    ! celerity at level 0
    if (allocated(dummy_Vector11)) deallocate(dummy_Vector11)
    allocate(dummy_Vector11(level0(L0_Basin(iBasin))%ncells))
    dummy_Vector11(:) = 0.0_dp
    call append(L0_celerity, dummy_Vector11)

    ! free space
    if (allocated(dummy_Vector11)) deallocate(dummy_Vector11)
    if (allocated(dummy_Matrix11_IT)) deallocate(dummy_Matrix11_IT)

  end subroutine variables_alloc_routing

END MODULE mo_mrm_init
