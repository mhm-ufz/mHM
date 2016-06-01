!> \file mo_mrm_init.f90

!> \brief Wrapper for initializing Routing.

!> \details Calling all routines to initialize all mRM variables

!> \authors Luis Samaniego, Rohini Kumar and Stephan Thober
!> \date Aug 2015

MODULE mo_mrm_init

  ! This module sets the river network characteristics and routing order.

  ! Written  Luis Samaniego, Mar 2005

  IMPLICIT NONE

  public :: mrm_init
  public :: variables_default_init_routing

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         mrm_init

  !     PURPOSE
  !>        \brief Initialize all mRM variables at all levels (i.e., L0, L1, and L11).
  !
  !>        \details Initialize all mRM variables at all levels (i.e., L0, L1, and L11)
  !>        either with default values or with values from restart file. The L0 mask (L0_mask),
  !>        L0 elevation (L0_elev), and L0 land cover (L0_LCover) can be provided as optional
  !>        variables to save memory because these variable will then not be read in again.
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !>        \param[in] "logical, dimension(:), target, optional :: L0_mask - L0 mask"
  !>        \param[in] "real(dp), dimension(:), target, optional :: L0_elev - L0 elevation"
  !>        \param[in] "integer(i4), dimension(:,:), target, optional :: L0_LCover - L0 land cover"
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !>        None
  !
  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Aug 2015
  !         Modified, Sep 2015 - Stephan Thober, added L0_mask, L0_elev, and L0_LCover
  !                   May 2016 - Stephan Thober, added warning message in case no gauge is found in modelling domain

  subroutine mrm_init(L0_mask, L0_elev, L0_LCover)

    use mo_kind,    only: i4, dp
    use mo_message, only: message
    use mo_mrm_constants, only: nodata_i4
    use mo_mrm_global_variables, only: read_restart, nBasins, perform_mpr, L0_Basin, dirRestartIn, &
         mrm_coupling_mode, basin_mrm
    use mo_mrm_net_startup, only: &
         L11_variable_init, &
         L11_flow_direction, &
         L11_set_network_topology, &
         L11_routing_order, &
         L11_link_location, &
         L11_set_drain_outlet_gauges, &
         L11_stream_features
    use mo_mrm_read_config, only: read_mrm_config_coupling, read_mrm_config
    use mo_mrm_read_data, only: mrm_read_discharge, mrm_read_L0_data, &
         mrm_L1_variable_init, &
         mrm_L0_variable_init, &
         mrm_read_total_runoff
    use mo_mrm_read_latlon, only: read_latlon
    use mo_mrm_restart, only: mrm_read_restart_config

    implicit none
    ! input variables
    logical, dimension(:), target, intent(in), optional :: L0_mask ! L0 mask
    real(dp), dimension(:), target, intent(in), optional :: L0_elev ! L0 elevation
    integer(i4), dimension(:,:), target, intent(in), optional :: L0_LCover ! L0 land cover

    ! local variables
    integer(i4) :: iBasin
    logical :: ReadLatLon

    !-----------------------------------------------------------
    ! READ COUPLING MODE
    !-----------------------------------------------------------
    call read_mrm_config_coupling()
    
    !-----------------------------------------------------------
    ! PRINT STARTUP MESSAGE
    !-----------------------------------------------------------
    if (mrm_coupling_mode .eq. 0) then
       call print_startup_message()
    else
       call message('')
       call message('  Inititalize mRM')
    end if
    
    ! ----------------------------------------------------------
    ! READ CONFIG
    ! ----------------------------------------------------------
    call read_mrm_config((mrm_coupling_mode .ne. 2), ReadLatLon)

    !-----------------------------------------------------------
    ! CONFIG OUTPUT
    !-----------------------------------------------------------
    if (mrm_coupling_mode .eq. 0) then
       call config_output()
    end if

    ! ----------------------------------------------------------
    ! READ DATA
    ! ----------------------------------------------------------
    ! level 0 data
    call mrm_read_L0_data(L0_mask, L0_elev, L0_LCover)
    
    if (perform_mpr) then
       do iBasin = 1, nBasins
          if (iBasin .eq. 1) then
             call L0_check_input_routing(iBasin)
          else if (L0_Basin(iBasin) .ne. L0_Basin(iBasin - 1) ) then
             call L0_check_input_routing(iBasin)
          end if
       end do
    end if

    if (read_restart) then
       ! ----------------------------------------------------------
       ! READ RESTART VARIABLES
       ! ----------------------------------------------------------
       do iBasin = 1, nBasins
          call mrm_read_restart_config(iBasin, dirRestartIn(iBasin))
       end do
    else
       ! ----------------------------------------------------------
       ! LEVEL 1 DATA
       ! ----------------------------------------------------------
       do iBasin = 1, nBasins
          if (iBasin .eq. 1) then
             call mrm_L0_variable_init(iBasin)
          else if (L0_Basin(iBasin) .ne. L0_Basin(iBasin - 1 )) then
             call mrm_L0_variable_init(iBasin)
          end if
          call mrm_L1_variable_init(iBasin)
       end do
       ! ----------------------------------------------------------
       ! INITIALIZE STREAM NETWORK
       ! ----------------------------------------------------------
       do iBasin = 1, nBasins
          call L11_variable_init(iBasin)
          call L11_flow_direction(iBasin)
          call L11_set_network_topology(iBasin)
          call L11_routing_order(iBasin)
          call L11_link_location(iBasin)
          call L11_set_drain_outlet_gauges(iBasin)
          ! stream characteristics
          call L11_stream_features(iBasin)
       end do
       ! check whether there are gauges within the modelling domain
       if (allocated(basin_mrm%gaugeNodeList)) then
          if (all(basin_mrm%gaugeNodeList .eq. nodata_i4)) then
             call message('')
             call message('    WARNING: no gauge found within modelling domain')
          end if
       end if
    end if

    ! ----------------------------------------------------------
    ! INITIALIZE STATES
    ! ----------------------------------------------------------
    do iBasin = 1, nBasins
       call variables_alloc_routing(iBasin)
    end do

    ! -------------------------------------------------------
    ! READ INPUT DATA AND OBSERVED DISCHARGE DATA
    ! -------------------------------------------------------
    ! read simulated runoff at level 1
    if (mrm_coupling_mode .eq. 0) then
       do iBasin = 1, nBasins
          call mrm_read_total_runoff(iBasin)
       end do
    end if
    ! discharge data
    call mrm_read_discharge()

    ! ----------------------------------------------------------
    ! READ LAT & LON IF OUTPUT SHALL BE WRITTEN
    ! ----------------------------------------------------------
    if (ReadLatLon) then
       do iBasin = 1, nBasins
          call read_latlon(iBasin)
       end do
    end if

    call message('')
    call message('  Finished Initialization of mRM')

  end subroutine mrm_init

  
  !===============================================================
  ! PRINT STARTUP MESSAGE
  !===============================================================
  subroutine print_startup_message()

    use mo_kind,    only: i4
    use mo_message, only: message, message_text
    use mo_string_utils, only: num2str, separator
    use mo_mrm_file, only: &
         version, &
         version_date, &
         file_main, &
         file_namelist_mrm, &
         file_namelist_param_mrm, &
         file_defOutput
    implicit none
    ! local variables
    integer(i4), dimension(8) :: datetime ! Date and time

    call message(separator)
    call message('              mRM-UFZ')
    call message()
    call message('    MULTISCALE ROUTING MODEL')
    call message('           Version ', trim(version))
    call message('           ', trim(version_date))
    call message()
    call message('Made available by S. Thober')
    call message()
    call message('Based on mHM-UFZ by L. Samaniego & R. Kumar')

    call message(separator)

    call message()
    call date_and_time(values=datetime)
    message_text = trim(num2str(datetime(3),'(I2.2)'))//"."//trim(num2str(datetime(2),'(I2.2)')) &
         //"."//trim(num2str(datetime(1),'(I4.4)'))//" "//trim(num2str(datetime(5),'(I2.2)')) &
         //":"//trim(num2str(datetime(6),'(I2.2)'))//":"//trim(num2str(datetime(7),'(I2.2)'))
    call message('Start at ', trim(message_text), '.')
    call message('Using main file ', trim(file_main), ' and namelists: ')
    call message('     ',trim(file_namelist_mrm))
    call message('     ',trim(file_namelist_param_mrm))
    call message('     ', trim(file_defOutput), ' (if it is given)')
    call message()

  end subroutine print_startup_message

  
  !---------------------------------------------------------------
  ! Config output
  !---------------------------------------------------------------
  subroutine config_output()

    use mo_kind,    only: i4
    use mo_string_utils, only: num2str
    use mo_message, only: message
    use mo_mrm_file, only: &
         file_namelist_mrm, &
         file_namelist_param_mrm, &
         file_defOutput
    use mo_mrm_global_variables, only: &
         dirMorpho, &
         dirLCover, &
         dirGauges, &
         dirOut, &
         basin_mrm, &
         nBasins
    implicit none
    ! local variables
    integer(i4) :: ii
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
    do ii = 1, nbasins
       call message( '  --------------' )
       call message( '      BASIN                   ', num2str(ii,'(I3)') )
       call message( '  --------------' )
       call message('    Morphological directory:    ',   trim(dirMorpho(ii) ))
       call message('    Land cover directory:       ',   trim(dirLCover(ii) ))
       call message('    Discharge directory:        ',   trim(dirGauges(ii)  ))
       call message('    Output directory:           ',   trim(dirOut(ii) ))
       call message('    Evaluation gauge            ', 'ID')
       do jj = 1 , basin_mrm%nGauges(ii)
          call message('    ',trim(adjustl(num2str(jj))),'                           ', &
               trim(adjustl(num2str(basin_mrm%gaugeIdList(ii,jj)))))
       end do
       if (basin_mrm%nInflowGauges(ii) .GT. 0) then
          call message('    Inflow gauge              ', 'ID')
          do jj = 1 , basin_mrm%nInflowGauges(ii)
             call message('    ',trim(adjustl(num2str(jj))),'                         ', &
                  trim(adjustl(num2str(basin_mrm%InflowGaugeIdList(ii,jj)))))
          end do
       end if
    end do
  end subroutine config_output

 
  ! ------------------------------------------------------------------

  !      NAME
  !          variables_default_init_routing

  !>        \brief Default initalization mRM related L11 variables
  
  !>        \details Default initalization of mHM related L11 variables (e.g., states,
  !>                 fluxes, and parameters) as per given constant values given in mo_mhm_constants.
  !>                 Variables initalized here is defined in the mo_global_variables.f90 file. 
  !>                 Only Variables that are defined in the variables_alloc subroutine are 
  !>                 intialized here.
  !
  !>                 If a variable is added or removed here, then it also has to be added or removed
  !>                 in the subroutine state_variables_set in the module mo_restart and in the 
  !>                 subroutine set_state in the module mo_set_netcdf_restart.
  !
  !     CALLING SEQUENCE
  !         call variables_default_init()

  !     INTENT(IN)

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL

  !     RETURN

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE

  !     HISTORY
  !         \authors  Stephan Thober, Rohini Kumar, and Juliane Mai
  !         \date    Aug 2015
  !         Modified

  subroutine variables_default_init_routing()

    use mo_mrm_global_variables, only: &
         L11_Qmod, L11_qOUT, L11_qTIN,  L11_qTR, L11_K, L11_xi,L11_C1, L11_C2,  &
         L11_FracFPimp
    use mo_mrm_constants,    only:               &
         P1_InitStateFluxes
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
    ! Fraction of the flood plain with impervious cover
    L11_FracFPimp = P1_InitStateFluxes
  end subroutine variables_default_init_routing

  
  ! --------------------------------------------------------------------------
  ! L0_check_input_routing
  ! --------------------------------------------------------------------------
  subroutine L0_check_input_routing(iBasin)
    
    use mo_kind,    only: i4
    use mo_message, only: message, message_text
    use mo_string_utils, only: num2str
    use mo_mrm_global_variables, only: L0_fDir, L0_fAcc, basin_mrm
    USE mo_mrm_constants, ONLY: nodata_i4
    implicit none
    ! input variables
    integer(i4) :: iBasin
    ! local variables
    integer(i4) :: k

    do k = basin_mrm%L0_iStart(iBasin), basin_mrm%L0_iEnd(iBasin)
       ! flow direction [-]
       if ( L0_fDir(k) .eq. nodata_i4  ) then
          message_text = trim(num2str(k,'(I5)'))//','// trim(num2str(iBasin,'(I5)'))
          call message(' Error: flow direction has missing value within the valid masked area at cell in basin ', &
               trim(message_text) )
          stop
       end if
       ! flow accumulation [-]
       if ( L0_fAcc(k) .eq. nodata_i4 ) then
          message_text = trim(num2str(k,'(I5)'))//','// trim(num2str(iBasin,'(I5)'))
          call message(' Error: flow accumulation has missing values within the valid masked area at cell in basin ', &
               trim(message_text) )
          stop
       end if
    end do
    
  end subroutine L0_check_input_routing

  
  ! --------------------------------------------------------------------------
  ! L11 ROUTING STATE VARIABLES, FLUXES AND
  !             PARAMETERS
  ! --------------------------------------------------------------------------
  subroutine variables_alloc_routing(iBasin)
    
    use mo_kind,    only: i4, dp
    use mo_mrm_constants,    only: nRoutingStates
    use mo_append,           only: append                      ! append vector
    use mo_mrm_tools, only: get_basin_info_mrm
    use mo_mrm_global_variables, only: L11_Qmod, L11_qOUT, L11_qTIN, &
         L11_qTR, L11_K, L11_xi,L11_C1, L11_C2, L11_FracFPimp
    implicit none
    ! input variables
    integer(i4), intent(in) :: iBasin
    ! local variables
    integer(i4) :: nrows11, ncols11
    integer(i4) :: ncells11
    real(dp), dimension(:),   allocatable :: dummy_Vector11
    real(dp), dimension(:,:), allocatable :: dummy_Matrix11_IT

    ! level-11 information
    call get_basin_info_mrm( iBasin, 11, nrows11, ncols11, ncells=ncells11 ) 

    ! dummy vector and matrix
    allocate( dummy_Vector11   (nCells11     ) )
    allocate( dummy_Matrix11_IT(nCells11, nRoutingStates) )

    ! simulated discharge at each node
    dummy_Vector11(:) = 0.0_dp
    call append( L11_Qmod,  dummy_Vector11 )

    ! Total outflow from cells L11 at time tt
    dummy_Vector11(:) = 0.0_dp
    call append( L11_qOUT, dummy_Vector11 )

    ! Total discharge inputs at t-1 and t
    dummy_Matrix11_IT(:,:) = 0.0_dp
    call append( L11_qTIN, dummy_Matrix11_IT )

    !  Routed outflow leaving a node
    dummy_Matrix11_IT(:,:) = 0.0_dp
    call append( L11_qTR, dummy_Matrix11_IT )

    ! kappa: Muskingum travel time parameter.
    dummy_Vector11(:) = 0.0_dp
    call append( L11_K, dummy_Vector11 )

    ! xi:    Muskingum diffusion parameter
    dummy_Vector11(:) = 0.0_dp
    call append( L11_xi, dummy_Vector11 )

    ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
    dummy_Vector11(:) = 0.0_dp
    call append( L11_C1, dummy_Vector11 )

    ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
    dummy_Vector11(:) = 0.0_dp
    call append( L11_C2, dummy_Vector11 )

    ! Fraction of the flood plain with impervious cover
    dummy_Vector11(:) = 0.0_dp
    call append( L11_FracFPimp, dummy_Vector11 )

    ! free space
    if ( allocated( dummy_Vector11        ) ) deallocate( dummy_Vector11        )
    if ( allocated( dummy_Matrix11_IT     ) ) deallocate( dummy_Matrix11_IT     )

  end subroutine variables_alloc_routing

END MODULE mo_mrm_init
