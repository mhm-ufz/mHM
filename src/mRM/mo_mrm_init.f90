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

  ! TODO: MPR no mrm_configuration, but mrm_update_param moved here
  public :: mrm_init, mrm_configuration

  private

CONTAINS

subroutine mrm_configuration(file_namelist, unamelist, file_namelist_param, unamelist_param)
    use mo_common_variables, only : processMatrix
    use mo_mrm_read_config, only : mrm_read_config
    use mo_mrm_global_variables, only: riv_temp_pcs
    use mo_common_read_config, only : check_optimization_settings
    use mo_kind, only : i4
    use mo_message, only : message
    implicit none

    character(*), intent(in) :: file_namelist, file_namelist_param

    integer, intent(in) :: unamelist, unamelist_param

    call message('')
    call message('  Inititalize mRM')
    if ( processMatrix(11, 1) .ne. 0 ) then
      ! processCase(11): river temperature routing
      riv_temp_pcs%active = .true.
      riv_temp_pcs%case = processMatrix(11, 1)
      call message('')
      call message('    Read config: river temperature routing')
      call riv_temp_pcs%config(file_namelist, unamelist, file_namelist_param, unamelist_param)
    end if

    ! read config for mrm, readlatlon is set here depending on whether output is needed
    call mrm_read_config(file_namelist, unamelist, file_namelist_param, unamelist_param, .false.)

end subroutine mrm_configuration
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
  ! Stephan Thober May 2019 - added init of level0 in case of read restart


  subroutine mrm_init(file_namelist, unamelist, file_namelist_param, unamelist_param)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_grid, only : read_grid_info
    use mo_common_variables, only : domainMeta, global_parameters, l0_l1_remap, level0, level1, domainMeta, &
                                    processMatrix, resolutionHydrology, mrmFileRestartIn, &
                                    mrm_read_river_network, resolutionRouting
    use mo_grid, only : init_lowres_level, set_domain_indices
    use mo_kind, only : i4
    use mo_message, only : message
    use mo_mrm_global_variables, only : domain_mrm, &
                                        l0_l11_remap, l1_l11_remap, level11, &
                                        gw_coupling, L0_river_head_mon_sum, &
                                        L11_netPerm, L11_fromN, L11_length, L11_nOutlets, &
                                        riv_temp_pcs, &
                                        readLatLon
    use mo_mrm_net_startup, only : L11_flow_direction, L11_flow_accumulation, L11_fraction_sealed_floodplain, &
                                   L11_link_location, L11_routing_order, L11_set_drain_outlet_gauges, &
                                   L11_set_network_topology, L11_stream_features, l11_l1_mapping
    use mo_mrm_read_data, only : mrm_read_L0_data, mrm_read_discharge, &
                                 mrm_read_total_runoff, mrm_read_bankfull_runoff
    use mo_mrm_restart, only : mrm_read_restart_config
    use mo_read_latlon, only : read_latlon
    use mo_mrm_river_head, only: init_masked_zeros_l0, create_output, calc_channel_elevation
    use mo_mrm_mpr, only : mrm_init_param
    use mo_common_read_config, only : check_optimization_settings

    implicit none

    character(*), intent(in) :: file_namelist, file_namelist_param

    integer, intent(in) :: unamelist, unamelist_param

    ! start and end index for routing parameters
    integer(i4) :: iStart, iEnd
    ! start and end index at L11
    integer(i4) :: s11, e11

    integer(i4) :: domainID, iDomain, gauge_counter

    call message('')
    call message('  Inititalize mRM')

    ! ----------------------------------------------------------
    ! READ DATA
    ! ----------------------------------------------------------
    allocate(level11(domainMeta%nDomains))
    allocate(l0_l1_remap(domainMeta%nDomains))
    allocate(l0_l11_remap(domainMeta%nDomains))
    allocate(l1_l11_remap(domainMeta%nDomains))

    if (.not. allocated(level0)) allocate(level0(domainMeta%nDomains))

    ! INIT all level 0 related information
    if (.not. mrm_read_river_network) then
      ! read all necessary level 0 data
      ! do_reinit, do_readlatlon, do_readlcover
      if (processMatrix(8, 1) .eq. 1_i4) call mrm_read_L0_data(.true.)
      if (processMatrix(8, 1) .eq. 2_i4) call mrm_read_L0_data(.false.)
      if (processMatrix(8, 1) .eq. 3_i4) call mrm_read_L0_data(.false.)
    else
      do iDomain = 1, domainMeta%nDomains
        domainID = domainMeta%indices(iDomain)
        ! this reads the domain properties
        ! ToDo: L0_Domain, parallel
        call read_grid_info(domainMeta%indices(domainMeta%L0DataFrom(iDomain)), mrmFileRestartIn(iDomain), &
                                                     "0", level0(domainMeta%L0DataFrom(iDomain)))
      end do
    end if
    call set_domain_indices(level0, indices=domainMeta%L0DataFrom)

    ! INIT all level 11 and 1 related information
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      if (mrm_read_river_network) then
        ! read information of l11 grids directly from restart file
        call read_grid_info(domainID, mrmFileRestartIn(iDomain), "11", level11(iDomain))
        call mrm_read_restart_config(iDomain, mrmFileRestartIn(iDomain))
      else
        ! TODO: this sets the lat lon values (can be projected!) and effective cellArea of level1
        ! think about how this can be done in MPR
        call init_lowres_level(level0(domainMeta%L0DataFrom(iDomain)), resolutionHydrology(iDomain), &
                level1(iDomain), l0_l1_remap(iDomain))
        ! init l11 grids from l0 grid
        call init_lowres_level(level0(domainMeta%L0DataFrom(iDomain)), resolutionRouting(iDomain), &
                level11(iDomain), l0_l11_remap(iDomain))
        call init_lowres_level(level1(iDomain), resolutionRouting(iDomain), &
                level11(iDomain), l1_l11_remap(iDomain))
        call L11_L1_mapping(iDomain)

        if (ReadLatLon) then
          ! read lat lon coordinates of each domain
          call read_latlon(iDomain, "lon_l11", "lat_l11", "level11", level11(iDomain))
        else
          ! allocate the memory and set to nodata
          level11(iDomain)%x = nodata_dp
          level11(iDomain)%y = nodata_dp
        end if
      end if
    end do

    call set_domain_indices(level11)
    call set_domain_indices(level1)

    ! ----------------------------------------------------------
    ! INITIALIZE STATES AND ROUTING PARAMETERS
    ! ----------------------------------------------------------
    do iDomain = 1, domainMeta%nDomains
      call variables_alloc_routing(iDomain)
    end do
    !-------------------------------------------
    ! L11 ROUTING STATE VARIABLES, FLUXES AND
    !             PARAMETERS
    !-------------------------------------------
    call variables_default_init_routing()

    ! ----------------------------------------------------------
    ! INITIALIZE STREAM NETWORK
    ! ----------------------------------------------------------
    do iDomain = 1, domainMeta%nDomains
       if (.not. mrm_read_river_network) then
        call L11_flow_direction(iDomain)
        call L11_flow_accumulation(iDomain)
        call L11_set_network_topology(iDomain)
        call L11_routing_order(iDomain)
        call L11_link_location(iDomain)
        call L11_set_drain_outlet_gauges(iDomain)
        ! stream characteristics
        call L11_stream_features(iDomain)
      end if

    end do

    ! check whether there are gauges within the modelling domain
    if (allocated(domain_mrm)) then
      gauge_counter = 0
      do iDomain = 1, domainMeta%nDomains
        if (.not. all(domain_mrm(iDomain)%gaugeNodeList .eq. nodata_i4)) then
          gauge_counter = gauge_counter + 1
        end if
      end do
      if (gauge_counter .lt. 1) then
        call message('')
        call message('    WARNING: no gauge found within modelling domain')
      end if
    end if
    ! ----------------------------------------------------------
    ! INITIALIZE PARAMETERS
    ! ----------------------------------------------------------
    do iDomain = 1, domainMeta%nDomains
      iStart = processMatrix(8, 3) - processMatrix(8, 2) + 1
      iEnd = processMatrix(8, 3)
      call mrm_init_param(iDomain, global_parameters(iStart : iEnd, 3))
    end do

    ! mpr-like definiton of sealed floodplain fraction
    if ((processMatrix(8, 1) .eq. 1_i4) .and. (.not. mrm_read_river_network)) then
      call L11_fraction_sealed_floodplain(2_i4, .true.)
    else
      ! dummy initialization
      call L11_fraction_sealed_floodplain(2_i4, .false.)
    end if

    ! -------------------------------------------------------
    ! READ INPUT DATA AND OBSERVED DISCHARGE DATA
    ! -------------------------------------------------------
    ! discharge data
    call mrm_read_discharge()

    if (gw_coupling) then
        do iDomain = 1, domainMeta%nDomains
            call init_masked_zeros_l0(iDomain, L0_river_head_mon_sum)
            call mrm_read_bankfull_runoff(iDomain)
            call create_output(iDomain, dirOut(iDomain))
        end do
        call calc_channel_elevation()
    end if

    ! init riv temp
    if ( riv_temp_pcs%active ) then
      call message('')
      call message('    Initialization of river temperature routing.')
      do iDomain = 1, domainMeta%nDomains
        s11 = level11(iDomain)%iStart
        e11 = level11(iDomain)%iEnd
        call riv_temp_pcs%init(level11(iDomain)%nCells)
        call riv_temp_pcs%init_area( &
          iDomain, &
          L11_netPerm(s11 : e11), & ! routing order at L11
          L11_fromN(s11 : e11), & ! link source at L11
          L11_length(s11 : e11 - 1), & ! link length
          level11(iDomain)%nCells - L11_nOutlets(iDomain), &
          level11(iDomain)%nCells, &
          level11(iDomain)%nrows, &
          level11(iDomain)%ncols, &
          level11(iDomain)%mask &
        )
      end do
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
    use mo_message, only : message
    use mo_mrm_file, only : file_defOutput, file_main, version, version_date
    use mo_string_utils, only : num2str, separator

    implicit none

    character(*), intent(in) :: file_namelist, file_namelist_param
    character(4096) :: message_text

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

    use mo_common_variables, only : dirOut, domainMeta
    use mo_kind, only : i4
    use mo_message, only : message
    use mo_mrm_file, only : file_defOutput, file_namelist_mrm, file_namelist_param_mrm
    use mo_mrm_global_variables, only : domain_mrm, &
                                        dirGauges
    use mo_string_utils, only : num2str

    implicit none

    integer(i4) :: domainID, iDomain

    integer(i4) :: jj


    !
    call message()
    call message('Read namelist file: ', trim(file_namelist_mrm))
    call message('Read namelist file: ', trim(file_namelist_param_mrm))
    call message('Read namelist file: ', trim(file_defOutput), ' (if it is given)')

    call message()
    call message('  # of domains:         ', trim(num2str(domainMeta%nDomains)))
    call message()
    call message('  Input data directories:')
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      call message('  --------------')
      call message('      DOMAIN                   ', num2str(domainID, '(I3)'))
      call message('  --------------')
      call message('    Discharge directory:        ', trim(dirGauges(iDomain)))
      call message('    Output directory:           ', trim(dirOut(iDomain)))
      call message('    Evaluation gauge            ', 'ID')
      do jj = 1, domain_mrm(iDomain)%nGauges
        call message('    ', trim(adjustl(num2str(jj))), '                           ', &
                trim(adjustl(num2str(domain_mrm(iDomain)%gaugeIdList(jj)))))
      end do
      if (domain_mrm(iDomain)%nInflowGauges .GT. 0) then
        call message('    Inflow gauge              ', 'ID')
        do jj = 1, domain_mrm(iDomain)%nInflowGauges
          call message('    ', trim(adjustl(num2str(jj))), '                         ', &
                  trim(adjustl(num2str(domain_mrm(iDomain)%InflowGaugeIdList(jj)))))
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
  ! L11 ROUTING STATE VARIABLES, FLUXES AND
  !             PARAMETERS
  ! --------------------------------------------------------------------------
  !    NAME
  !        variables_alloc_routing

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:

  subroutine variables_alloc_routing(iDomain)

    use mo_append, only : append
    use mo_kind, only : dp, i4
    use mo_mrm_constants, only : nRoutingStates
    use mo_common_variables, only : level0, domainMeta
    use mo_mrm_global_variables, only : L11_C1, L11_C2, L11_K, &
         L11_Qmod, L11_qOUT, L11_qTIN, L11_qTR, L11_xi, &
         level11, L11_celerity, L0_celerity

    implicit none

    integer(i4), intent(in) :: iDomain

    real(dp), dimension(:), allocatable :: dummy_Vector11

    real(dp), dimension(:, :), allocatable :: dummy_Matrix11_IT


    ! dummy vector and matrix
    allocate(dummy_Vector11   (level11(iDomain)%nCells))
    allocate(dummy_Matrix11_IT(level11(iDomain)%nCells, nRoutingStates))

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
    allocate(dummy_Vector11(level0(domainMeta%L0DataFrom(iDomain))%ncells))
    dummy_Vector11(:) = 0.0_dp
    call append(L0_celerity, dummy_Vector11)

    ! free space
    if (allocated(dummy_Vector11)) deallocate(dummy_Vector11)
    if (allocated(dummy_Matrix11_IT)) deallocate(dummy_Matrix11_IT)

  end subroutine variables_alloc_routing

END MODULE mo_mrm_init
