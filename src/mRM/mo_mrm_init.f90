!> \file mo_mrm_init.f90
!> \brief \copybrief mo_mrm_init
!> \details \copydetails mo_mrm_init

!> \brief Wrapper for initializing Routing.
!> \details Calling all routines to initialize all mRM variables
!> \authors Luis Samaniego, Rohini Kumar and Stephan Thober
!> \date Aug 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
MODULE mo_mrm_init

    use mo_common_variables, only : dirOut
    use mo_message, only : message, error_message

  ! This module sets the river network characteristics and routing order.

  ! Written  Luis Samaniego, Mar 2005

  IMPLICIT NONE

  public :: mrm_init, mrm_configuration
  public :: variables_default_init_routing
  public :: fluxes_states_default_init_routing

  private

CONTAINS

  !> \brief read mRM configuration from namelists
  subroutine mrm_configuration(file_namelist, unamelist, file_namelist_param, unamelist_param)
    use mo_common_mHM_mRM_variables, only : mrm_coupling_mode
    use mo_common_variables, only : processMatrix
    use mo_mrm_read_config, only : mrm_read_config
    use mo_mrm_global_variables, only: riv_temp_pcs
    use mo_common_read_config, only : common_read_config
    use mo_common_mHM_mRM_read_config, only : check_optimization_settings, common_mHM_mRM_read_config
    use mo_kind, only : i4
    implicit none

    character(*), intent(in) :: file_namelist !< namelist file name
    integer, intent(in) :: unamelist !< unit to open namelist
    character(*), intent(in) :: file_namelist_param !< parameter namelist file name
    integer, intent(in) :: unamelist_param !< unit to open parameter namelist

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
      if ( processMatrix(11, 1) .ne. 0 ) then
        ! processCase(11): river temperature routing
        riv_temp_pcs%active = .true.
        riv_temp_pcs%case = processMatrix(11, 1)
        call message('')
        call message('    Read config: river temperature routing')
        call riv_temp_pcs%config(file_namelist, unamelist, file_namelist_param, unamelist_param)
      end if
    end if

    ! read config for mrm, readlatlon is set here depending on whether output is needed
    call mrm_read_config(file_namelist, unamelist, file_namelist_param, unamelist_param, (mrm_coupling_mode .eq. 0_i4))

    ! this was moved here, because it depends on global_parameters that are only set in mrm_read_config
    if (mrm_coupling_mode .eq. 0_i4) then
      call check_optimization_settings()
      !-----------------------------------------------------------
      ! CONFIG OUTPUT
      !-----------------------------------------------------------
      call config_output()
    end if
  end subroutine mrm_configuration


  !> \brief Initialize all mRM variables at all levels (i.e., L0, L1, and L11).
  !> \details Initialize all mRM variables at all levels (i.e., L0, L1, and L11)
  !! either with default values or with values from restart file. The L0 mask (L0_mask),
  !! L0 elevation (L0_elev), and L0 land cover (L0_LCover) can be provided as optional
  !! variables to save memory because these variable will then not be read in again.
  !> \changelog
  !! - Stephan Thober Sep 2015
  !!   - added L0_mask, L0_elev, and L0_LCover
  !! - Stephan Thober May 2016
  !!   - added warning message in case no gauge is found in modelling domain
  !! - Matthias Kelbling Aug 2017
  !!   - added L11_flow_accumulation to Initialize Stream Netwo
  !! - Lennart Schueler May 2018
  !!   - added initialization for groundwater coupling
  !! - Stephan Thober Jun 2018
  !!   - refactored for mpr_extract version
  !! - Stephan Thober May 2019
  !!   - added init of level0 in case of read restart
  !> \authors Stephan Thober
  !> \date Aug 2015
  subroutine mrm_init(file_namelist, unamelist, file_namelist_param, unamelist_param)

    use mo_common_constants, only : nodata_dp, nodata_i4
    use mo_common_mHM_mRM_variables, only : mrmFileRestartIn, mrm_coupling_mode, mrm_read_river_network, &
                                            resolutionRouting
    use mo_common_restart, only : read_grid_info
    use mo_common_variables, only : global_parameters, l0_l1_remap, level0, level1, domainMeta, &
                                    processMatrix, resolutionHydrology
    use mo_grid, only : L0_grid_setup, init_lowres_level, set_domain_indices
    use mo_kind, only : i4
    use mo_mrm_global_variables, only : domain_mrm, &
                                        l0_l11_remap, level11, &
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
    use mo_mrm_river_head, only: init_masked_zeros_l0, calc_channel_elevation
    use mo_mrm_mpr, only : mrm_init_param

    implicit none

    character(*), intent(in) :: file_namelist !< namelist file name
    integer, intent(in) :: unamelist !< unit to open namelist
    character(*), intent(in) :: file_namelist_param !< parameter namelist file name
    integer, intent(in) :: unamelist_param !< unit to open parameter namelist

    ! start and end index for routing parameters
    integer(i4) :: iStart, iEnd
    ! start and end index at L11
    integer(i4) :: s11, e11

    integer(i4) :: domainID, iDomain, gauge_counter


    if (mrm_coupling_mode .eq. 0_i4) then
      allocate(l0_l1_remap(domainMeta%nDomains))
      allocate(level1(domainMeta%nDomains))
    end if

    ! ----------------------------------------------------------
    ! READ DATA
    ! ----------------------------------------------------------
    allocate(level11(domainMeta%nDomains))
    allocate(l0_l11_remap(domainMeta%nDomains))

    if (.not. mrm_read_river_network) then
      ! read all (still) necessary level 0 data
      if (processMatrix(8, 1) .eq. 1_i4) call mrm_read_L0_data(mrm_coupling_mode .eq. 0_i4, ReadLatLon, .true.)
      if (processMatrix(8, 1) .eq. 2_i4) call mrm_read_L0_data(mrm_coupling_mode .eq. 0_i4, ReadLatLon, .false.)
      if (processMatrix(8, 1) .eq. 3_i4) call mrm_read_L0_data(mrm_coupling_mode .eq. 0_i4, ReadLatLon, .false.)
    end if

    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      if (mrm_read_river_network) then
        ! this reads the domain properties
        if (.not. allocated(level0)) allocate(level0(domainMeta%nDomains))
        ! ToDo: L0_Domain, parallel
        call read_grid_info(mrmFileRestartIn(iDomain), "0", level0(domainMeta%L0DataFrom(iDomain)))
        if (mrm_coupling_mode .eq. 0_i4) then
          call read_grid_info(mrmFileRestartIn(iDomain), "1", level1(iDomain))
        end if
        call read_grid_info(mrmFileRestartIn(iDomain), "11", level11(iDomain))
        call mrm_read_restart_config(iDomain, domainID, mrmFileRestartIn(iDomain))
      else
        if (iDomain .eq. 1) then
          call L0_check_input_routing(domainMeta%L0DataFrom(iDomain))
          if (mrm_coupling_mode .eq. 0_i4) then
            call L0_grid_setup(level0(domainMeta%L0DataFrom(iDomain)))
          end if
        else if ((domainMeta%L0DataFrom(iDomain) == iDomain)) then
          call L0_check_input_routing(domainMeta%L0DataFrom(iDomain))
          if (mrm_coupling_mode .eq. 0_i4) then
            call L0_grid_setup(level0(domainMeta%L0DataFrom(iDomain)))
          end if
        end if

        if (mrm_coupling_mode .eq. 0_i4) then
          call init_lowres_level(level0(domainMeta%L0DataFrom(iDomain)), resolutionHydrology(iDomain), &
                  level1(iDomain), l0_l1_remap(iDomain))
        end if
        call init_lowres_level(level0(domainMeta%L0DataFrom(iDomain)), resolutionRouting(iDomain), &
                level11(iDomain), l0_l11_remap(iDomain))
        call L11_L1_mapping(iDomain)

        if (ReadLatLon) then
          ! read lat lon coordinates of each domain
          call read_latlon(iDomain, "lon", "lat", "level1", level1(iDomain))
          call read_latlon(iDomain, "lon_l11", "lat_l11", "level11", level11(iDomain))
        else
          ! allocate the memory and set to nodata
          allocate(level11(iDomain)%x(level11(iDomain)%nrows, level11(iDomain)%ncols))
          allocate(level11(iDomain)%y(level11(iDomain)%nrows, level11(iDomain)%ncols))
          level11(iDomain)%x = nodata_dp
          level11(iDomain)%y = nodata_dp
        end if
      end if
    end do

    call set_domain_indices(level11)
    call set_domain_indices(level1)
    call set_domain_indices(level0, indices=domainMeta%L0DataFrom)

    ! ----------------------------------------------------------
    ! INITIALIZE STATES AND AUXILLIARY VARIABLES
    ! ----------------------------------------------------------
    do iDomain = 1, domainMeta%nDomains
      call variables_alloc_routing(iDomain)
    end do

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

    ! ----------------------------------------------------------
    ! INITIALIZE PARAMETERS
    ! ----------------------------------------------------------
    do iDomain = 1, domainMeta%nDomains
      iStart = processMatrix(8, 3) - processMatrix(8, 2) + 1
      iEnd = processMatrix(8, 3)
      call mrm_init_param(iDomain, global_parameters(iStart : iEnd, 3))
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
    ! read simulated runoff at level 1
    if (mrm_coupling_mode .eq. 0_i4) then
      do iDomain = 1, domainMeta%nDomains
        call mrm_read_total_runoff(iDomain)
      end do
    end if
    ! discharge data
    call mrm_read_discharge()

    ! init groundwater coupling
    if (gw_coupling) then
      do iDomain = 1, domainMeta%nDomains
        call init_masked_zeros_l0(iDomain, L0_river_head_mon_sum)
        call mrm_read_bankfull_runoff(iDomain)
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


  !> \brief Print mRM startup message
  !> \authors Robert Schweppe
  !> \date Jun 2018
  subroutine print_startup_message(file_namelist, file_namelist_param)

    use mo_kind, only : i4
    use mo_mrm_file, only : file_defOutput, file_main, version, version_date
    use mo_string_utils, only : num2str, separator

    implicit none

    character(*), intent(in) :: file_namelist !< namelist file name
    character(*), intent(in) :: file_namelist_param !< parameter namelist file name

    ! Date and time
    integer(i4), dimension(8) :: datetime

    CHARACTER(len=1024) :: message_text = ''

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


  !> \brief print mRM configuration
  !> \authors Robert Schweppe
  !> \date Jun 2018
  subroutine config_output

    use mo_common_variables, only : dirLCover, dirMorpho, dirOut, domainMeta
    use mo_kind, only : i4
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
      call message('    Morphological directory:    ', trim(dirMorpho(iDomain)))
      call message('    Land cover directory:       ', trim(dirLCover(iDomain)))
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


  !> \brief Default initalization mRM related L11 variables
  !> \details Default initalization of mHM related L11 variables (e.g., states,
  !! fluxes, and parameters) as per given constant values given in mo_mhm_constants.
  !! Variables initalized here is defined in the mo_global_variables.f90 file.
  !! Only Variables that are defined in the variables_alloc subroutine are
  !! intialized here.
  !! If a variable is added or removed here, then it also has to be added or removed
  !! in the subroutine state_variables_set in the module mo_restart and in the
  !! subroutine set_state in the module mo_set_netcdf_restart.
  !> \authors  Stephan Thober, Rohini Kumar, and Juliane Mai
  !> \date    Aug 2015
  !> \authors Robert Schweppe
  !> \date Jun 2018
  subroutine variables_default_init_routing

    use mo_common_constants, only : P1_InitStateFluxes
    use mo_mrm_global_variables, only : L11_C1, L11_C2, L11_K, L11_xi

    implicit none

    !-------------------------------------------
    ! L11 ROUTING STATE VARIABLES, FLUXES AND
    !             PARAMETERS
    !-------------------------------------------

    ! fluxes and states
    call fluxes_states_default_init_routing()

    ! kappa: Muskingum travel time parameter.
    L11_K = P1_InitStateFluxes
    ! xi:    Muskingum diffusion parameter
    L11_xi = P1_InitStateFluxes
    ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
    L11_C1 = P1_InitStateFluxes
    ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
    L11_C2 = P1_InitStateFluxes

  end subroutine variables_default_init_routing

  !> \brief initialize fluxes and states with default values for mRM
  subroutine fluxes_states_default_init_routing(iDomain)

    use mo_kind, only: i4
    use mo_mrm_global_variables, only : level11
    use mo_common_constants, only : P1_InitStateFluxes
    use mo_mrm_global_variables, only : L11_Qmod, L11_qOUT, L11_qTIN, L11_qTR

    implicit none

    !> number of Domain (if not present, set for all)
    integer(i4), intent(in), optional :: iDomain

    integer(i4) :: s11, e11

    !-------------------------------------------
    ! L11 ROUTING STATE VARIABLES, FLUXES AND
    !             PARAMETERS
    !-------------------------------------------

    if (present(iDomain)) then
      s11 = level11(iDomain)%iStart
      e11 = level11(iDomain)%iEnd
      ! simulated discharge at each node
      L11_Qmod(s11 : e11) = P1_InitStateFluxes
      ! Total outflow from cells L11 at time tt
      L11_qOUT(s11 : e11) = P1_InitStateFluxes
      ! Total discharge inputs at t-1 and t
      L11_qTIN(s11 : e11, :) = P1_InitStateFluxes
      !  Routed outflow leaving a node
      L11_qTR(s11 : e11, :) = P1_InitStateFluxes
    else
      ! simulated discharge at each node
      L11_Qmod = P1_InitStateFluxes
      ! Total outflow from cells L11 at time tt
      L11_qOUT = P1_InitStateFluxes
      ! Total discharge inputs at t-1 and t
      L11_qTIN = P1_InitStateFluxes
      !  Routed outflow leaving a node
      L11_qTR = P1_InitStateFluxes
    end if

  end subroutine fluxes_states_default_init_routing


  !> \brief check routing input on level-0
  !> \authors Robert Schweppe
  !> \date Jun 2018
  subroutine L0_check_input_routing(L0Domain_iDomain)

    use mo_common_constants, only : nodata_i4
    use mo_common_variables, only : level0
    use mo_kind, only : i4
    use mo_mrm_global_variables, only : L0_fAcc, L0_fDir
    use mo_string_utils, only : num2str

    implicit none

    integer(i4), intent(in) :: L0Domain_iDomain !< domain index for associated level-0 data

    integer(i4) :: k

    CHARACTER(len=1024) :: message_text = ''

    do k = level0(L0Domain_iDomain)%iStart, level0(L0Domain_iDomain)%iEnd
      ! flow direction [-]
      if (L0_fDir(k) .eq. nodata_i4) then
        message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(L0Domain_iDomain, '(I5)'))
        call error_message(' Error: flow direction has missing value within the valid masked area at cell in domain ', &
                trim(message_text))
      end if
      ! flow accumulation [-]
      if (L0_fAcc(k) .eq. nodata_i4) then
        message_text = trim(num2str(k, '(I5)')) // ',' // trim(num2str(L0Domain_iDomain, '(I5)'))
        call error_message(' Error: flow accumulation has missing values within the valid masked area at cell in domain ', &
                trim(message_text))
      end if
    end do

  end subroutine L0_check_input_routing


  !> \brief allocated routing related variables
  !> \authors Robert Schweppe
  !> \date Jun 2018
  subroutine variables_alloc_routing(iDomain)

    use mo_append, only : append
    use mo_kind, only : dp, i4
    use mo_mrm_constants, only : nRoutingStates
    use mo_common_variables, only : level0, domainMeta
    use mo_mrm_global_variables, only : L11_C1, L11_C2, L11_K, &
         L11_Qmod, L11_qOUT, L11_qTIN, L11_qTR, L11_xi, &
         level11, L11_celerity, L0_celerity

    implicit none

    integer(i4), intent(in) :: iDomain !< domain index

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
