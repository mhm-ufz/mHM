!> \file mo_restart.f90

!> \brief reading and writing states, fluxes and configuration for restart of mHM.

!> \details routines are seperated for reading and writing variables for:\n
!>          - states and fluxes, and \n
!>          - configuration.\n
!>          Reading of L11 configuration is also seperated from the rest, 
!>          since it is only required when routing is activated.

!> \authors Stephan Thober
!> \date Jul 2013

MODULE mo_restart

  ! This module is a restart for the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Stephan Thober, Apr 2011

  IMPLICIT NONE

  PUBLIC :: read_restart_states     ! read restart files for state variables from a given path
  PUBLIC :: read_restart_config     ! read restart files for configuration from a given path
  PUBLIC :: read_restart_L11_config ! read L11 configuration
  PUBLIC :: write_restart_files     ! write restart files for configuration to a given path

  PRIVATE

CONTAINS
  ! ------------------------------------------------------------------
  
  !      NAME
  !         write_restart

  !     PURPOSE
  !>        \brief write restart files for each basin

  !>        \details write restart files for each basin. For each basin
  !>        three restart files are written. These are xxx_states.nc, 
  !>        xxx_L11_config.nc, and xxx_config.nc (xxx being the three digit
  !>        basin index). If a variable is added here, it should also be added
  !>        in the read restart routines below.

  !     INTENT(IN)
  !>        \param[in] "character(256), dimension(:) :: OutPath"     Output Path for each basin

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN

  !     RESTRICTIONS 
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         see library routine var2nc in mo_ncwrite.f90

  !     HISTORY
  !>        \author   Stephan Thober
  !>        \date     Jun 2014
  !         Modified   Matthias Zink       Nov. 2014  - added PET related parameter writing
  ! ------------------------------------------------------------------ 
  subroutine write_restart_files( OutPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_init_states,      only: get_basin_info
    use mo_string_utils,     only: num2str
    use mo_ncwrite,          only: var2nc
    use mo_mhm_constants,    only: nodata_dp, nodata_i4, nRoutingStates
    use mo_global_variables, only: processMatrix, &
         L1_fSealed, &
         L1_fForest, &
         L1_fPerm, &
         L1_Inter, &
         L1_snowPack, &
         L1_sealSTW, &
         L1_soilMoist, &
         L1_unsatSTW, &
         L1_satSTW, &
         L1_aETSoil, &
         L1_aETCanopy, &
         L1_aETSealed, &
         L1_baseflow, &
         L1_infilSoil, &
         L1_fastRunoff, &
         L1_melt, &
         L1_percol, &
         L1_preEffect, &
         L1_rain, &
         L1_runoffSeal, &
         L1_slowRunoff, &
         L1_snow, &
         L1_Throughfall, &
         L1_total_runoff, &
         L1_alpha, &
         L1_degDayInc, &
         L1_degDayMax, &
         L1_degDayNoPre, &
         L1_degDay, &
         L1_karstLoss, &
         L1_fAsp, &
         L1_HarSamCoeff, &
         L1_PrieTayAlpha, &
         L1_aeroResist, &
         L1_surfResist, &
         L1_fRoots, &
         L1_maxInter, &
         L1_kfastFlow, &
         L1_kSlowFlow, &
         L1_kBaseFlow, &
         L1_kPerco, &
         L1_soilMoistFC, &
         L1_soilMoistSat, &
         L1_soilMoistExp, &
         L1_tempThresh, &
         L1_unsatThresh, &
         L1_sealedThresh, &
         L1_wiltingPoint, &
         L11_Qmod, &
         L11_qOUT, &
         L11_qTIN, &
         L11_qTR, &
         L11_K, &
         L11_xi, &
         L11_C1, &
         L11_C2, &
         L11_FracFPimp, &
         L11_cellCoor, &
         L11_Id, &
         L0_L11_Id, &
         L1_L11_Id, &
         L11_rowOut, &
         L11_colOut, &
         L11_fDir, &
         L11_upBound_L0, &
         L11_downBound_L0, &
         L11_leftBound_L0, &
         L11_rightBound_L0, &
         L11_upBound_L1, &
         L11_downBound_L1, &
         L11_leftBound_L1, &
         L11_rightBound_L1, &
         L11_fDir, &
         L11_fromN, &
         L11_toN, &
         L11_rOrder, &
         L11_label, &
         L11_sink, &
         L11_netPerm, &
         L11_rowOut, &
         L11_colOut, &
         L11_fRow, &
         L11_fCol, &
         L11_tRow, &
         L11_tCol, &
         L0_draSC, &
         L0_draCell, &
         L0_streamNet, &
         L0_floodPlain, &
         L11_length, &
         L11_aFloodPlain, &
         L11_slope, &
         basin, &
         L0_cellCoor    ,          & 
         L0_Id         ,           & ! Ids of grid at level-0 
         L0_areaCell   ,           & ! Ids of grid at level-0
         L0_slope_emp  ,           & ! Empirical quantiles of slope
         L1_Id         ,           & ! Ids of grid at level-1
         L1_cellCoor    ,          &
         L1_upBound_L0 ,           & ! Row start at finer level-0 scale 
         L1_downBound_L0,          & ! Row end at finer level-0 scale
         L1_leftBound_L0,          & ! Col start at finer level-0 scale
         L1_rightBound_L0,         & ! Col end at finer level-0 scale
         L1_areaCell   ,           & ! [km2] Effective area of cell at this level
         L1_nTCells_L0               ! Total number of valid L0 cells in a given L1 cell

    implicit none

    character(256)                           :: Fname
    character(256), dimension(:), intent(in) :: OutPath ! list of Output paths per Basin
    integer(i4)                              :: iBasin
    integer(i4)                              :: ii
    integer(i4)                              :: s0       ! start index at level 0
    integer(i4)                              :: e0       ! end index at level 0
    integer(i4)                              :: ncols0   ! number of colums at level 0
    integer(i4)                              :: nrows0   ! number of rows at level 0
    logical, dimension(:,:), allocatable     :: mask0    ! mask at level 0
    integer(i4)                              :: s1       ! start index at level 1
    integer(i4)                              :: e1       ! end index at level 1
    integer(i4)                              :: ncols1   ! number of colums at level 1
    integer(i4)                              :: nrows1   ! number of rows at level 1
    logical, dimension(:,:), allocatable     :: mask1    ! mask at level 1
    integer(i4)                              :: s110     ! start index at pseudo level 110 
    integer(i4)                              :: e110     ! end index at pseudo level 110
    integer(i4)                              :: ncols110 ! number of colums at pseudo level 110
    integer(i4)                              :: nrows110 ! number of rows at pseudo level 110
    integer(i4)                              :: s11      ! start index at level 11
    integer(i4)                              :: e11      ! end index at level 11
    integer(i4)                              :: ncols11  ! number of colums at level 11
    integer(i4)                              :: nrows11  ! number of rows at level 11
    logical, dimension(:,:), allocatable     :: mask11   ! mask at level 11
    real(dp), dimension(:,:,:), allocatable  :: dummy_d3 ! dummy variable
    ! dimension variables
    character(256), dimension(2)             :: dims_L0     ! dimension names for L0 states
    character(256), dimension(4)             :: dims_L1     ! dimension names for L1 states
    character(256), dimension(3)             :: dims_L11    ! dimension names for L11 states
    character(256), dimension(1)             :: dims_outlet ! dimension name  for outlet Coordinates
    character(256), dimension(1)             :: dims_gauges ! dimension name  for number of gauges
    character(256), dimension(1)             :: dims_inflow ! dimension name  for inflow gauge

    ! initialize
    dims_L0(1)     = 'nrows0'
    dims_L0(2)     = 'ncols0'
    dims_L1(1)     = 'nrows1'
    dims_L1(2)     = 'ncols1'
    dims_L1(3)     = 'L1_soilhorizons'
    dims_L1(4)     = 'MonthsPerYear'
    dims_L11(1)    = 'nrows11'
    dims_L11(2)    = 'ncols11'
    dims_L11(3)    = 'nIT'
    dims_outlet(1) = 'NoutletCoord'
    dims_gauges(1) = 'Ngauges'
    dims_inflow(1) = 'nInflowGauges'

    basin_loop: do iBasin = 1, size(OutPath)

       ! get Level0 information about the basin
       call get_basin_info( iBasin, 0, nrows0, ncols0, iStart=s0, iEnd=e0, mask=mask0 )

       ! get Level1 information about the basin
       call get_basin_info( iBasin, 1, nrows1, ncols1, iStart=s1, iEnd=e1, mask=mask1 )

       ! write restart file for iBasin
       ! Fname = trim(OutPath(iBasin)) // trim(num2str(iBasin, '(i3.3)')) // '_restart.nc'
       Fname = trim(OutPath(iBasin)) // trim(num2str(iBasin, '(i3.3)')) // '_states.nc'
       ! print a message
       call message('    Writing Restart-file: ', trim(adjustl(Fname)),' ...')
       
       call var2nc( Fname, unpack( L1_fSealed(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_fSealed', &
            long_name = 'fraction of Sealed area at level 1', missing_value = nodata_dp, &
            create = .true. ) ! create file

       call var2nc( Fname, unpack( L1_fForest(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_fForest', &
            long_name = 'fraction of Forest area at level 1', missing_value = nodata_dp)
       
       call var2nc( Fname, unpack( L1_fPerm(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_fPerm', &
            long_name = 'fraction of permeable area at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_inter(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_Inter', &
            long_name = 'Interception storage at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_snowPack(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_snowPack', &
            long_name = 'Snowpack at level 1', missing_value = nodata_dp)
       
       call var2nc( Fname, unpack( L1_sealSTW(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_sealSTW', &
            long_name = 'Retention storage of impervious areas at level 1', missing_value = nodata_dp)

       allocate( dummy_d3( nrows1, ncols1, size( L1_soilMoist, 2) ) )
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_soilMoist(s1:e1,ii), mask1, nodata_dp )
       end do
       call var2nc( Fname, dummy_d3, &
            dims_L1(1:3), 'L1_soilMoist', &
            long_name = 'soil moisture at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_unsatSTW(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_unsatSTW', &
            long_name = 'upper soil storage at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_satSTW(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_satSTW', &
            long_name = 'groundwater storage at level 1', missing_value = nodata_dp)

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_aETSoil(s1:e1,ii), mask1, nodata_dp )
       end do
       call var2nc( Fname, dummy_d3, &
            dims_L1(1:3), 'L1_aETSoil', &
            long_name = 'soil actual ET at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_aETCanopy(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_aETCanopy', &
            long_name = 'canopy actual ET at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_aETSealed(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_aETSealed', &
            long_name = 'sealed actual ET at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_baseflow(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_baseflow', &
            long_name = 'baseflow at level 1', missing_value = nodata_dp)

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_infilSoil(s1:e1,ii), mask1, nodata_dp )
       end do
       call var2nc( Fname, dummy_d3, &
            dims_L1(1:3), 'L1_infilSoil', &
            long_name = 'soil in-exfiltration at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_fastRunoff(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_fastRunoff', &
            long_name = 'fast runoff', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_percol(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_percol', &
            long_name = 'percolation at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_melt(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_melt', &
            long_name = 'snow melt at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_preEffect(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_preEffect', &
            long_name = 'effective precip. depth (snow melt + rain) at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_rain(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_rain', &
            long_name = 'rain (liquid water) at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_runoffSeal(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_runoffSeal', &
            long_name = 'runoff from impervious area at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_slowRunoff(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_slowRunoff', &
            long_name = 'slow runoff at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_snow(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_snow', &
            long_name = 'snow (solid water) at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_Throughfall(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_Throughfall', &
            long_name = 'throughfall at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_total_runoff(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_total_runoff', &
            long_name = 'total runoff at level 1', missing_value = nodata_dp)

       !-------------------------------------------
       ! EFFECTIVE PARAMETERS
       !-------------------------------------------
       call var2nc( Fname, unpack( L1_alpha(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_alpha', &
            long_name = 'exponent for the upper reservoir at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_degDayInc(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_degDayInc', &
            long_name = 'increase of the Degree-day factor per mm of increase in precipitation at level 1', &
            missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_degDayMax(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_degDayMax', &
            long_name = 'maximum degree-day factor at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_degDayNoPre(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_degDayNoPre', &
            long_name = 'degree-day factor with no precipitation at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_degDay(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_degDay', &
            long_name = 'degree-day factor at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_karstLoss(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_karstLoss', &
            long_name = 'Karstic percolation loss at level 1', missing_value = nodata_dp)

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_fRoots(s1:e1,ii), mask1, nodata_dp )
       end do
       call var2nc( Fname, dummy_d3, &
            dims_L1(1:3), 'L1_fRoots', &
            long_name = 'Fraction of roots in soil horizons at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_maxInter(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_maxInter', &
            long_name = 'Maximum interception at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_kfastFlow(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_kfastFlow', &
            long_name = 'fast interflow recession coefficient at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_kSlowFlow(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_kSlowFlow', &
            long_name = 'slow interflow recession coefficient at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_kBaseFlow(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_kBaseFlow', &
            long_name = 'baseflow recession coefficient at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_kPerco(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_kPerco', &
            long_name = 'percolation coefficient at level 1', missing_value = nodata_dp)

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_soilMoistFC(s1:e1,ii), mask1, nodata_dp )
       end do
       call var2nc( Fname, dummy_d3, &
            dims_L1(1:3), 'L1_soilMoistFC', &
            long_name = 'Soil moisture below which actual ET is reduced linearly till PWP at level 1', missing_value = nodata_dp)

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_soilMoistSat(s1:e1,ii), mask1, nodata_dp )
       end do
       call var2nc( Fname, dummy_d3, &
            dims_L1(1:3), 'L1_soilMoistSat', &
            long_name = 'Saturation soil moisture for each horizon [mm] at level 1', missing_value = nodata_dp)

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_soilMoistExp(s1:e1,ii), mask1, nodata_dp )
       end do
       call var2nc( Fname, dummy_d3, &
            dims_L1(1:3), 'L1_soilMoistExp', &
            long_name = 'Exponential parameter to how non-linear is the soil water retention at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_tempThresh(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_tempThresh', &
            long_name = 'Threshold temperature for snow/rain at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_unsatThresh(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_unsatThresh', &
            long_name = 'Threshhold water depth controlling fast interflow at level 1', missing_value = nodata_dp)

       call var2nc( Fname, unpack( L1_sealedThresh(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2), 'L1_sealedThresh', &
            long_name = 'Threshhold water depth for surface runoff in sealed surfaces at level 1', missing_value = nodata_dp)

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_wiltingPoint(s1:e1,ii), mask1, nodata_dp )
       end do
       call var2nc( Fname, dummy_d3, &
            dims_L1(1:3), 'L1_wiltingPoint', &
            long_name = 'Permanent wilting point at level 1', missing_value = nodata_dp)
       
       deallocate( dummy_d3 )

       select case (processMatrix(5,1))
       case(0) ! PET is input
          call var2nc( Fname, unpack( L1_fAsp(s1:e1), mask1, nodata_dp ), &
               dims_L1(1:2), 'L1_fAsp', &
               long_name = 'PET correction factor due to terrain aspect at level 1', missing_value = nodata_dp)        
       case(1) ! HarSam
          call var2nc( Fname, unpack( L1_fAsp(s1:e1), mask1, nodata_dp ), &
               dims_L1(1:2), 'L1_fAsp', &
               long_name = 'PET correction factor due to terrain aspect at level 1', missing_value = nodata_dp)        
          call var2nc( Fname, unpack( L1_HarSamCoeff(s1:e1), mask1, nodata_dp ), &
               dims_L1(1:2), 'L1_HarSamCoeff', &
               long_name = 'Hargreaves-Samani coefficient', missing_value = nodata_dp)        
       case(2) ! PrieTay
          allocate( dummy_d3( nrows1, ncols1, size( L1_PrieTayAlpha, 2) ) )
          do ii = 1, size( dummy_d3, 3 )
             dummy_d3(:,:,ii) = unpack( L1_PrieTayAlpha(s1:e1,ii), mask1, nodata_dp )
          end do

          call var2nc( Fname, dummy_d3, &
               (/dims_L1(1:2),dims_L1(4)/), 'L1_PrieTayAlpha', &
               long_name = 'Priestley Taylor coeffiecient (alpha)', missing_value = nodata_dp)        
          deallocate( dummy_d3 )
       case(3) ! PenMon
          allocate( dummy_d3( nrows1, ncols1, size( L1_aeroResist, 2) ) )
          do ii = 1, size( dummy_d3, 3 )
             dummy_d3(:,:,ii) = unpack( L1_aeroResist(s1:e1,ii), mask1, nodata_dp )
          end do
          call var2nc( Fname, dummy_d3, &
               (/dims_L1(1:2),dims_L1(4)/), 'L1_aeroResist', &
               long_name = 'aerodynamical resitance', missing_value = nodata_dp)        

          do ii = 1, size( dummy_d3, 3 )
             dummy_d3(:,:,ii) = unpack( L1_surfResist(s1:e1,ii), mask1, nodata_dp )
          end do
          call var2nc( Fname, dummy_d3, &
               (/dims_L1(1:2),dims_L1(4)/), 'L1_surfResist', &
               long_name = 'bulk surface resitance', missing_value = nodata_dp)        
          deallocate( dummy_d3 )
       end select

       !-------------------------------------------
       ! L11 ROUTING STATE VARIABLES, FLUXES AND
       !             PARAMETERS
       !-------------------------------------------
       if ( processMatrix(8,1) .ne. 0 ) then
          ! get Level11 information about the basin
          call get_basin_info( iBasin, 11, nrows11, ncols11, iStart=s11, iEnd=e11, mask=mask11 )
          ! get Level110 information about the basin
          call get_basin_info( iBasin, 110, nrows110, ncols110, iStart=s110, iEnd=e110)
          
          call var2nc( Fname, unpack( L11_Qmod(s11:e11), mask11, nodata_dp ), &
            dims_L11(1:2), 'L11_Qmod', &
            long_name = 'simulated discharge at each node at level 11', missing_value = nodata_dp)

          call var2nc( Fname, unpack( L11_qOUT(s11:e11), mask11, nodata_dp ), &
            dims_L11(1:2), 'L11_qOUT', &
            long_name = 'Total outflow from cells L11 at time tt at level 11', missing_value = nodata_dp)
          
          allocate( dummy_d3( nrows11, ncols11, nRoutingStates ) )
          do ii = 1, size( dummy_d3, 3 )
             dummy_d3(:,:,ii) = unpack( L11_qTIN(s11:e11,ii), mask11, nodata_dp )
          end do
          call var2nc( Fname, dummy_d3, &
               dims_L11, 'L11_qTIN', &
               long_name = 'Total discharge inputs at t-1 and t at level 11', missing_value = nodata_dp)

          do ii = 1, size( dummy_d3, 3 )
             dummy_d3(:,:,ii) = unpack( L11_qTR(s11:e11,ii), mask11, nodata_dp )
          end do
          call var2nc( Fname, dummy_d3, &
               dims_L11, 'L11_qTR', &
               long_name = 'Routed outflow leaving a node at level 11', missing_value = nodata_dp)

          call var2nc( Fname, unpack( L11_K(s11:e11), mask11, nodata_dp ), &
            dims_L11(1:2), 'L11_K', &
            long_name = 'kappa: Muskingum travel time parameter at level 11', missing_value = nodata_dp)

          call var2nc( Fname, unpack( L11_xi(s11:e11), mask11, nodata_dp ), &
            dims_L11(1:2), 'L11_xi', &
            long_name = 'xi: Muskingum diffusion parameter at level 11', missing_value = nodata_dp)

          call var2nc( Fname, unpack( L11_C1(s11:e11), mask11, nodata_dp ), &
            dims_L11(1:2), 'L11_C1', &
            long_name = 'Routing parameter C1=f(K,xi, DT) (Chow, 25-41) at level 11', missing_value = nodata_dp)

          call var2nc( Fname, unpack( L11_C2(s11:e11), mask11, nodata_dp ), &
            dims_L11(1:2), 'L11_C2', &
            long_name = 'Routing parameter C2=f(K,xi, DT) (Chow, 25-41) at level 11', missing_value = nodata_dp)

          call var2nc( Fname, unpack( L11_FracFPimp(s11:e11), mask11, nodata_dp ), &
            dims_L11(1:2), 'L11_FracFPimp', &
            long_name = 'Fraction of the flood plain with impervious cover at level 11', missing_value = nodata_dp)
          
          ! ----------------------------------------------------------
          ! L11 config set - create new file
          ! ----------------------------------------------------------
          Fname = trim(OutPath(iBasin)) // trim(num2str(iBasin, '(i3.3)')) // '_L11_config.nc'
          call message('    Writing Restart-file: ', trim(adjustl(Fname)),' ...')
          call var2nc( Fname, &
               merge( 1_i4, 0_i4,  &
               reshape(basin%L11_Mask(basin%L11_iStartMask(iBasin):basin%L11_iEndMask(iBasin)),&
               (/nrows11,ncols11/)) ),&
               dims_L11(1:2), 'L11_basin_Mask', &
               long_name = 'Mask at Level 11', missing_value = nodata_i4, create = .true. )
          
          call var2nc( Fname, unpack( L11_cellCoor(s11:e11,1), mask11, nodata_i4 ), &
               dims_L11(1:2),'L11_rowCoor', &
               long_name = 'row coordinates at Level 11', missing_value = nodata_i4 )

          call var2nc( Fname, unpack( L11_cellCoor(s11:e11,2), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_colCoor', &
               long_name = 'col coordinates at Level 11', missing_value = nodata_i4 )

          call var2nc( Fname, unpack( L11_Id(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_Id', &
               long_name = 'cell Ids at Level 11', missing_value = nodata_i4 )

          call var2nc( Fname, unpack( L11_fDir(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_fDir', &
               long_name = 'flow Direction at Level 11', missing_value = nodata_i4 )     

          call var2nc( Fname, unpack( L11_rowOut(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_rowOut', &
               long_name = 'Grid vertical location of the Outlet at Level 11', missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_colOut(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_colOut', &
               long_name = 'Grid horizontal location of the Outlet at Level 11',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_upBound_L0(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_upBound_L0', &
               long_name = 'Row start at finer level-0 scale of Level 11 cell',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_downBound_L0(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_downBound_L0', &
               long_name = 'Row end at finer level-0 scale of Level 11 cell',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_leftBound_L0(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_leftBound_L0', &
               long_name = 'Col start at finer level-0 scale of Level 11 cell',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_rightBound_L0(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_rightBound_L0', &
               long_name = 'Col end at finer level-0 scale of Level 11 cell',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_fromN(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_fromN', &
               long_name = 'From Node',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_toN(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_toN', &
               long_name = 'To Node',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_rOrder(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_rOrder', &
               long_name = 'Network routing order at Level 11',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_label(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_label', &
               long_name = 'Label Id [0='', 1=HeadWater, 2=Sink] at Level 11',missing_value=nodata_i4)

          call var2nc( Fname, unpack( merge( 1_i4, 0_i4, L11_sink(s11:e11)), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_sink', &
               long_name = '.true. if sink node reached at Level 11',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_netPerm(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_netPerm', &
               long_name = 'Routing sequence (permutation of L11_rOrder) at Level 11',missing_value=nodata_i4)
      
          call var2nc( Fname, unpack( L11_fRow(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_fRow', &
               long_name = 'From row in L0 grid at Level 11',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_fCol(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_fCol', &
               long_name = 'From col in L0 grid at Level 11',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_tRow(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_tRow', &
               long_name = 'To row in L0 grid at Level 11',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_tCol(s11:e11), mask11, nodata_i4 ), &
               dims_L11(1:2), 'L11_tCol', &
               long_name = 'To Col in L0 grid at Level 11',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_length(s11:e11), mask11, nodata_dp ), &
               dims_L11(1:2), 'L11_length', &
               long_name = 'Total length of river link [m]',missing_value=nodata_dp)

          call var2nc( Fname, unpack( L11_aFloodPlain(s11:e11), mask11, nodata_dp ), &
               dims_L11(1:2), 'L11_aFloodPlain', &
               long_name = 'Area of the flood plain [m2]',missing_value=nodata_dp)

          call var2nc( Fname, unpack( L11_slope(s11:e11), mask11, nodata_dp ), &
               dims_L11(1:2), 'L11_slope', &
               long_name = 'Average slope of river link',missing_value=nodata_dp)

          call var2nc( Fname, unpack( L0_draCell(s110:e110), mask0, nodata_i4 ), &
               dims_L0, 'L0_draCell', &
               long_name = 'Draining cell id at L11 of ith cell of L0',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L0_streamNet(s110:e110), mask0, nodata_i4 ), &
               dims_L0, 'L0_streamNet', &
               long_name = 'Stream network',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L0_floodPlain(s110:e110), mask0, nodata_i4 ), &
               dims_L0, 'L0_floodPlain', &
               long_name = 'Floodplains of stream i',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L0_draSC(s110:e110), mask0, nodata_i4 ), &
               dims_L0, 'L0_draSC', &
               long_name = 'Floodplains of stream i',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L0_L11_Id(s110:e110), mask0, nodata_i4 ), &
               dims_L0, 'L0_L11_Id', &
               long_name = 'Mapping of L11 Id on L0',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L1_L11_Id(s1:e1), mask1, nodata_i4 ), &
               dims_L1(1:2), 'L1_L11_Id', &
               long_name = 'Mapping of L11 Id on L1',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_upBound_L1(s1:e1), mask1, nodata_i4 ), &
               dims_L1(1:2), 'L11_upBound_L1', &
               long_name = 'Row start at finer level-1 scale',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_downBound_L1(s1:e1), mask1, nodata_i4 ), &
               dims_L1(1:2), 'L11_downBound_L1', &
               long_name = 'Row end at finer level-1 scale',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_leftBound_L1(s1:e1), mask1, nodata_i4 ), &
               dims_L1(1:2), 'L11_leftBound_L1', &
               long_name = 'Col start at finer level-1 scale',missing_value=nodata_i4)

          call var2nc( Fname, unpack( L11_rightBound_L1(s1:e1), mask1, nodata_i4 ), &
               dims_L1(1:2), 'L11_rightBound_L1', &
               long_name = 'Col start at finer level-1 scale',missing_value=nodata_i4)

          call var2nc( Fname, (/ basin%L0_rowOutlet(iBasin), basin%L0_colOutlet(iBasin) /), &
                dims_outlet, 'L0_OutletCoord', &
                long_name = 'Outlet Coordinates at Level 0',missing_value=nodata_i4)
          
          call var2nc( Fname, basin%gaugeNodeList(iBasin,:), &
               dims_gauges, 'gaugeNodeList', &
               long_name = 'cell ID of gauges',missing_value=nodata_i4)

          call var2nc( Fname, basin%InflowGaugeNodeList(iBasin,:), &
               dims_inflow, 'InflowGaugeNodeList', &
               long_name = 'cell ID of gauges',missing_value=nodata_i4)
       
          ! free dummy variables
          deallocate( dummy_d3 )

       end if

       ! -------------------------------------------------------------
       ! config set - create new file
       ! -------------------------------------------------------------
       Fname = trim(OutPath(iBasin)) // trim(num2str(iBasin, '(i3.3)')) // '_config.nc'
       call message('    Writing Restart-file: ', trim(adjustl(Fname)),' ...')

       call var2nc( Fname, unpack( L0_cellCoor(s0:e0,1), mask0, nodata_i4 ), &
            dims_L0,'L0_rowCoor', &
            long_name = 'row coordinates at Level 0', missing_value = nodata_i4, &
            create = .true. )

       call var2nc( Fname, unpack( L0_cellCoor(s0:e0,2), mask0, nodata_i4 ), &
            dims_L0,'L0_colCoor', &
            long_name = 'col coordinates at Level 0', missing_value = nodata_i4 )

       call var2nc( Fname, unpack( L0_Id(s0:e0), mask0, nodata_i4 ), &
            dims_L0,'L0_Id', &
            long_name = 'cell IDs at level 0', missing_value = nodata_i4 )

       call var2nc( Fname, unpack( L0_areaCell(s0:e0), mask0, nodata_dp ), &
            dims_L0,'L0_areaCell', &
            long_name = 'Area of a cell at level-0 [m2]', missing_value = nodata_dp )

       call var2nc( Fname, unpack( L0_slope_emp(s0:e0), mask0, nodata_dp ), &
            dims_L0,'L0_slope_emp', &
            long_name = 'Empirical quantiles of slope', missing_value = nodata_dp )
 
       call var2nc( Fname, &
            merge( 1_i4, 0_i4,  &
            reshape(basin%L1_Mask(basin%L1_iStartMask(iBasin):basin%L1_iEndMask(iBasin)),&
            (/nrows1,ncols1/)) ),&
            dims_L1(1:2), 'L1_basin_Mask', &
            long_name = 'Mask at Level 1', missing_value = nodata_i4 )
 
       call var2nc( Fname, unpack( L1_Id(s1:e1), mask1, nodata_i4 ), &
            dims_L1(1:2),'L1_Id', &
            long_name = 'cell IDs at level 1', missing_value = nodata_i4 )

       call var2nc( Fname, unpack( L1_cellCoor(s1:e1,1), mask1, nodata_i4 ), &
            dims_L1(1:2),'L1_rowCoor', &
            long_name = 'row cell Coordinates at Level 1', missing_value = nodata_i4 )

       call var2nc( Fname, unpack( L1_cellCoor(s1:e1,2), mask1, nodata_i4 ), &
            dims_L1(1:2),'L1_colCoor', &
            long_name = 'col cell Coordinates at Level 1', missing_value = nodata_i4 )
 
       call var2nc( Fname, unpack( L1_upBound_L0(s1:e1), mask1, nodata_i4 ), &
            dims_L1(1:2),'L1_upBound_L0', &
            long_name = 'Row start at finer level-0 scale', missing_value = nodata_i4 )
 
       call var2nc( Fname, unpack( L1_downBound_L0(s1:e1), mask1, nodata_i4 ), &
            dims_L1(1:2),'L1_downBound_L0', &
            long_name = 'Row end at finer level-0 scale', missing_value = nodata_i4 )
  
       call var2nc( Fname, unpack( L1_leftBound_L0(s1:e1), mask1, nodata_i4 ), &
            dims_L1(1:2),'L1_leftBound_L0', &
            long_name = 'Col start at finer level-0 scale', missing_value = nodata_i4 )
 
       call var2nc( Fname, unpack( L1_rightBound_L0(s1:e1), mask1, nodata_i4 ), &
            dims_L1(1:2),'L1_rightBound_L0', &
            long_name = 'Col end at finer level-0 scal', missing_value = nodata_i4 )
 
       call var2nc( Fname, unpack( L1_nTCells_L0(s1:e1), mask1, nodata_i4 ), &
            dims_L1(1:2),'L1_nTCells_L0', &
            long_name = 'Total number of valid L0 cells in a given L1 cell', missing_value = nodata_i4 )
  
       call var2nc( Fname, unpack( L1_areaCell(s1:e1), mask1, nodata_dp ), &
            dims_L1(1:2),'L1_areaCell', &
            long_name = 'Effective area of cell at this level [km2]', missing_value = nodata_dp )
       
    end do basin_loop
    
    
  end subroutine write_restart_files
  ! ------------------------------------------------------------------

  !      NAMEw
  !         read_restart_L11_config

  !     PURPOSE
  !>        \brief reads Level 11 configuration from a restart directory

  !>        \details read Level 11 configuration variables from a given restart
  !>        directory and initializes all Level 11 configuration variables,
  !>        that are initialized in the subroutine initialise,
  !>        contained in module mo_startup.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"        number of basin
  !>        \param[in] "character(256) :: InPath"        Input Path including trailing slash

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN

  !     RESTRICTIONS 
  !>        \note Restart Files must have the format, as if
  !>        it would have been written by subroutine write_restart_files 

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author   Stephan Thober
  !>        \date     Apr 2013

  !         Modified  Matthias Zink , Apr 2014 - added inflow gauge

  subroutine read_restart_L11_config( iBasin, InPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_kind,             only: i4, dp
    use mo_init_states,      only: get_basin_info
    use mo_append,           only: append
    use mo_ncread,           only: Get_NcVar
    use mo_mhm_constants,    only: nodata_dp
    use mo_init_states,      only: calculate_grid_properties
    use mo_global_variables, only: basin, & ! basin database
         level11,           &
         resolutionRouting, &
         nBasins,           & ! Number of Basins
         L11_cellCoor,      & ! cell Coordinates at Level 11
         L11_Id,            & ! cell Ids at Level 11
         L11_nCells,        & ! Number of Cells at Level 11
         L0_draSC,          &
         L0_L11_Id,         &
         L1_L11_Id,         &
         L11_fDir,          &
         L11_rowOut,        &
         L11_colOut,        &
         L11_upBound_L0,    &
         L11_downBound_L0,  &
         L11_leftBound_L0,  &
         L11_rightBound_L0, &
         L11_upBound_L1,    &
         L11_downBound_L1,  &
         L11_leftBound_L1,  &
         L11_rightBound_L1, &
         L11_fromN,         &
         L11_toN,           &
         L11_rOrder,        &
         L11_label,         &
         L11_sink,          &
         L11_netPerm,       &
         L11_fRow,          &
         L11_fCol,          &
         L11_tRow,          &
         L11_tCol,          &
         L0_draCell,        &
         L0_streamNet,      &
         L0_floodPlain,     &
         L11_length,        &
         L11_aFloodPlain,   &
         L11_slope

    implicit none

    integer(i4),    intent(in) :: iBasin
    character(256), intent(in) :: InPath ! list of Output paths per Basin

    ! local variables
    character(256)             :: Fname

    ! local variables
    integer(i4)                                          :: nrows0   ! Number of rows at level 0
    integer(i4)                                          :: ncols0   ! Number of cols at level 0
    logical, dimension(:,:), allocatable                 :: mask0    ! Mask at Level 0
    integer(i4)                                          :: nrows1   ! Number of rows at level 1
    integer(i4)                                          :: ncols1   ! Number of cols at level 1
    logical, dimension(:,:), allocatable                 :: mask1    ! Mask at Level 1
    integer(i4)                                          :: nrows11  ! Number of rows at level 11
    integer(i4)                                          :: ncols11  ! Number of cols at level 11
    logical, dimension(:,:), allocatable                 :: mask11   ! Mask at Level 11
    integer(i4)                                          :: nCells11 ! Number of cells at lev. 11
    real(dp)                                             :: xllcorner0, yllcorner0
    real(dp)                                             :: cellsize0

    ! DUMMY variables
    integer(i4), dimension(:),     allocatable           :: dummyI1  ! dummy, 1 dimension I4
    integer(i4), dimension(:,:),   allocatable           :: dummyI2  ! dummy, 2 dimension I4
    integer(i4), dimension(:,:),   allocatable           :: dummyI22 ! 2nd dummy, 2 dimension I4
    real(dp),    dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension DP

    ! set file name
    Fname = trim(InPath) // trim(num2str(iBasin, '(i3.3)')) // '_L11_config.nc' ! '_restart.nc'
    call message('    Reading L11_config from ', trim(adjustl(Fname)),' ...')

    ! level-0 information
    call get_basin_info( iBasin, 0, nrows0, ncols0,&
         xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0, mask=mask0 )

    ! level-1 information
    call get_basin_info( iBasin, 1, nrows1, ncols1, mask=mask1 )

    ! calculate l11 grid resolutionRouting
    if(iBasin .eq. 1) then
       allocate( level11%nrows        (nBasins) )
       allocate( level11%ncols        (nBasins) )
       allocate( level11%xllcorner    (nBasins) )
       allocate( level11%yllcorner    (nBasins) )
       allocate( level11%cellsize     (nBasins) )
       allocate( level11%nodata_value (nBasins) )
    end if
    call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp,            &
         resolutionRouting(iBasin) , &
         level11%nrows(iBasin), level11%ncols(iBasin), level11%xllcorner(iBasin), &
         level11%yllcorner(iBasin), level11%cellsize(iBasin), level11%nodata_value(iBasin)        )

    ! level-11 information
    call get_basin_info (iBasin, 11, nrows11, ncols11)

    ! read L11 mask
    allocate( dummyI2( nrows11, ncols11 ), mask11( nrows11, ncols11) )
    call Get_NcVar( Fname,  'L11_basin_Mask', dummyI2 )
    mask11 = (dummyI2 .eq. 1_i4)
    call append( basin%L11_Mask, reshape( mask11, (/nrows11*ncols11/)))

    ! get Number of cells
    nCells11 = count( dummyI2 .eq. 1_i4 )
    deallocate( dummyI2 )

    ! update basin database
    if (iBasin .eq. 1) then

       !
       allocate(basin%L11_iStart     (nBasins))
       allocate(basin%L11_iEnd       (nBasins))
       allocate(basin%L11_iStartMask (nBasins))
       allocate(basin%L11_iEndMask   (nBasins))    

       ! basin information
       basin%L11_iStart(iBasin) = 1
       basin%L11_iEnd  (iBasin) = basin%L11_iStart(iBasin) + nCells11 - 1

       basin%L11_iStartMask(iBasin) = 1
       basin%L11_iEndMask  (iBasin) = basin%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    else

       ! basin information
       basin%L11_iStart(iBasin) = basin%L11_iEnd(iBasin-1) + 1
       basin%L11_iEnd  (iBasin) = basin%L11_iStart(iBasin) + nCells11 - 1

       basin%L11_iStartMask(iBasin) = basin%L11_iEndMask(iBasin-1) + 1
       basin%L11_iEndMask  (iBasin) = basin%L11_iStartMask(iBasin) + nrows11*ncols11 - 1

    end if

    ! read L11 cellCoor
    allocate( dummyI2( nrows11, ncols11 ) )
    allocate( dummyI22( count(mask11), 2 ))
    call Get_NcVar( Fname, 'L11_rowCoor', dummyI2 )
    dummyI22(:,1) = pack( dummyI2, mask11 )
    call Get_NcVar( Fname, 'L11_colCoor', dummyI2 )
    dummyI22(:,2) = pack( dummyI2, mask11 )
    call append( L11_cellCoor, dummyI22)
    deallocate( dummyI22 )

    ! read L11 IDs
    call Get_NcVar( Fname, 'L11_Id', dummyI2 )
    call append( L11_Id, pack( dummyI2, mask11) )
    deallocate( dummyI2 )

    ! update Number of cells at Level 11
    L11_nCells = size( L11_Id, 1 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L0 outlet Coordinates
    allocate(dummyI1(2))
    call Get_NcVar( Fname, 'L0_OutletCoord', dummyI1)

    ! allocate space for row and col Outlet
    if (iBasin .eq. 1) then
       allocate( basin%L0_rowOutlet(nBasins) ) 
       allocate( basin%L0_colOutlet(nBasins) )
    end if

    ! L0 data sets
    basin%L0_rowOutlet(iBasin) = dummyI1(1)
    basin%L0_colOutlet(iBasin) = dummyI1(2)
    deallocate(dummyI1)

    ! read L0 draining cell index
    allocate(dummyI2(nrows0,ncols0))
    call Get_NcVar( Fname, 'L0_draSC', dummyI2)
    call append( L0_draSC,     PACK ( dummyI2,  mask0)  )

    ! Mapping of L11 Id on L0
    call Get_NcVar( Fname, 'L0_L11_Id', dummyI2 )
    call append( L0_L11_Id,    PACK ( dummyI2, mask0)  )
    deallocate( dummyI2 )

    ! L1 data sets
    ! Mapping of L11 Id on L1
    allocate( dummyI2( nrows1, ncols1 ) )
    call Get_NcVar( Fname, 'L1_L11_Id', dummyI2 )
    call append( L1_L11_Id,    PACK ( dummyI2, mask1)  )
    deallocate( dummyI2 )

    ! L11 data sets
    ! Flow direction (standard notation)
    allocate( dummyI2( nrows11, ncols11) )
    call Get_NcVar( Fname, 'L11_fDir', dummyI2 )
    call append( L11_fDir,     PACK ( dummyI2,      mask11) )
    deallocate( dummyI2 )

    ! Grid vertical location of the Outlet
    allocate( dummyI2( nrows11, ncols11 ) )
    call Get_NcVar( Fname, 'L11_rowOut', dummyI2 )
    call append( L11_rowOut, pack( dummyI2, mask11 ) )

    ! Grid horizontal location  of the Outlet
    call Get_NcVar( Fname, 'L11_colOut', dummyI2 )
    call append( L11_colOut, pack( dummyI2, mask11 ) )

    ! Row start at finer level-0 scale
    call Get_NcVar( Fname, 'L11_upBound_L0', dummyI2 )
    call append( L11_upBound_L0, pack( dummyI2, mask11 ) )

    ! Row end at finer level-0 scale
    call Get_NcVar( Fname, 'L11_downBound_L0', dummyI2 )
    call append( L11_downBound_L0, pack( dummyI2, mask11) )

    ! Col start at finer level-0 scale
    call Get_NcVar( Fname, 'L11_leftBound_L0', dummyI2 )
    call append( L11_leftBound_L0, pack( dummyI2, mask11) )

    ! Col end at finer level-0 scale 
    call Get_NcVar( Fname, 'L11_rightBound_L0', dummyI2 )
    call append( L11_rightBound_L0, pack( dummyI2, mask11) )
    deallocate( dummyI2 )

    ! Row start at finer level-1 scale
    allocate( dummyI2( nrows1, ncols1 ) )
    call Get_NcVar( Fname, 'L11_upBound_L1', dummyI2 )
    call append( L11_upBound_L1, pack( dummyI2, mask1) )

    ! Row end at finer level-1 scale
    call Get_NcVar( Fname, 'L11_downBound_L1', dummyI2 )
    call append( L11_downBound_L1, pack( dummyI2, mask1) )

    ! Col start at finer level-1 scale
    call Get_NcVar( Fname, 'L11_leftBound_L1', dummyI2 )
    call append( L11_leftBound_L1, pack( dummyI2, mask1) )

    ! Col end at finer level-1 scale 
    call Get_NcVar( Fname, 'L11_rightBound_L1', dummyI2 )
    call append( L11_rightBound_L1, pack(dummyI2, mask1) ) 
    deallocate( dummyI2 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! From Node
    allocate( dummyI2( nrows11, ncols11 ) )
    call Get_NcVar( Fname, 'L11_fromN', dummyI2)
    call append( L11_fromN, pack( dummyI2, mask11) )

    ! To Node
    call Get_NcVar( Fname, 'L11_toN', dummyI2 )
    call append( L11_toN, pack( dummyI2, mask11) )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Network routing order
    call Get_NcVar( Fname, 'L11_rOrder', dummyI2 )
    call append( L11_rOrder, pack(dummyI2, mask11) )

    ! Label Id [0='', 1=HeadWater, 2=Sink]
    call Get_NcVar( Fname, 'L11_label', dummyI2 )
    call append( L11_label, pack( dummyI2, mask11) )

    ! .true. if sink node reached
    call Get_NcVar( Fname, 'L11_sink', dummyI2 )
    call append( L11_sink, (pack(dummyI2,mask11) .eq. 1_i4) )

    ! Routing sequence (permutation of L11_rOrder)
    call Get_NcVar( Fname, 'L11_netPerm', dummyI2 )
    call append( L11_netPerm, pack( dummyI2, mask11) )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! From row in L0 grid
    call Get_NcVar( Fname, 'L11_fRow', dummyI2 )
    call append( L11_fRow, pack( dummyI2, mask11) )

    ! From col in L0 grid
    call Get_NcVar( Fname, 'L11_fCol', dummyI2 )
    call append( L11_fCol, pack( dummyI2, mask11) )

    ! To row in L0 grid
    call Get_NcVar( Fname, 'L11_tRow', dummyI2 )
    call append( L11_tRow, pack( dummyI2, mask11) )

    ! To col in L0 grid
    call Get_NcVar( Fname, 'L11_tCol', dummyI2 )
    call append( L11_tCol, pack( dummyI2, mask11) )
    deallocate( dummyI2 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    allocate( dummyI2( nrows0, ncols0 ) )
    call Get_NcVar( Fname, 'L0_draCell', dummyI2 )
    call append( L0_draCell,     PACK ( dummyI2,  mask0)  ) 

    ! read gaugenodelist
    allocate(dummyI1( size(basin%gaugeNodeList(iBasin,:))))
    call Get_NcVar( Fname, 'gaugeNodeList', dummyI1)
    basin%gaugeNodeList( iBasin, : ) = dummyI1
    deallocate(dummyI1)

    ! read InflowGaugeNodelist
    if (basin%nInflowGauges(iBasin) > 0) then 
       allocate(dummyI1( size(basin%InflowGaugeNodeList(iBasin,:))))
       call Get_NcVar( Fname, 'InflowGaugeNodeList', dummyI1)
       basin%InflowgaugeNodeList( iBasin, : ) = dummyI1
       deallocate(dummyI1)
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! L0 data sets 
    ! Stream network
    call Get_NcVar( Fname, 'L0_streamNet', dummyI2 )
    call append( L0_streamNet,    PACK ( dummyI2,  mask0)   ) 

    ! Floodplains of stream i
    call Get_NcVar( Fname, 'L0_floodPlain', dummyI2 )
    call append( L0_floodPlain,   PACK ( dummyI2,  mask0)  ) 
    deallocate( dummyI2 )

    ! L11 network data sets
    ! [m]     Total length of river link
    allocate( dummyD2( nrows11, ncols11 ) )
    call Get_NcVar( Fname, 'L11_length', dummyD2 )
    call append( L11_length, pack( dummyD2, mask11 ) )

    ! [m2]    Area of the flood plain
    call Get_NcVar( Fname, 'L11_aFloodPlain', dummyD2 )
    call append( L11_aFloodPlain, pack( dummyD2, mask11) )

    ! Average slope of river link
    call Get_NcVar( Fname, 'L11_slope', dummyD2 )
    call append( L11_slope, pack( dummyD2, mask11 ) )
    deallocate( dummyD2 )

  end subroutine read_restart_L11_config

  ! ------------------------------------------------------------------

  !      NAME
  !         read_restart_config

  !     PURPOSE
  !>        \brief reads configuration apart from Level 11 configuration
  !>        from a restart directory

  !>        \details read configuration variables from a given restart
  !>        directory and initializes all configuration variables,
  !>        that are initialized in the subroutine initialise,
  !>        contained in module mo_startup.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"        number of basin
  !>        \param[in] "character(256) :: InPath"        Input Path including trailing slash

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN

  !     RESTRICTIONS 
  !>        \note Restart Files must have the format, as if
  !>        it would have been written by subroutine write_restart_files

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Apr 2013

  subroutine read_restart_config( iBasin, soilId_isPresent, InPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_append,           only: append
    use mo_ncread,           only: Get_NcVar
    use mo_mhm_constants,    only: nodata_dp
    use mo_init_states,      only: calculate_grid_properties
    use mo_global_variables, only: L0_Basin, & ! check whether L0_Basin should be read
         perform_mpr,    & ! switch that controls whether mpr is performed or not
         L0_soilId,      & ! soil IDs at lower level
         L0_cellCoor   , & 
         L0_Id         , & ! Ids of grid at level-0 
         L0_areaCell   , & ! Ids of grid at level-0
         L0_slope_emp  , & ! Empirical quantiles of slope
         basin, & 
         nBasins, &
         level1, &
         L0_nCells, &
         nSoilTypes, &
         resolutionHydrology, &
         L1_nCells,      &
         L1_Id         , & ! Ids of grid at level-1
         L1_cellCoor   , &
         L1_upBound_L0 , & ! Row start at finer level-0 scale 
         L1_downBound_L0, & ! Row end at finer level-0 scale
         L1_leftBound_L0, & ! Col start at finer level-0 scale
         L1_rightBound_L0, & ! Col end at finer level-0 scale
         L1_areaCell   , & ! [km2] Effective area of cell at this level
         L1_nTCells_L0     ! Total number of valid L0 cells in a given L1 cell

    implicit none

    !
    integer(i4), intent(in) :: iBasin
    integer(i4), dimension(:), allocatable, intent(inout) :: soilId_isPresent
    character(256), intent(in)     :: InPath             ! list of Output paths per Basin

    !
    integer(i4)                                          :: nrows0   ! Number of rows at level 0
    integer(i4)                                          :: ncols0   ! Number of cols at level 
    integer(i4)                                          :: iStart0, iEnd0
    logical, dimension(:,:), allocatable                 :: mask0    ! Mask at Level 0
    integer(i4)                                          :: nrows1   ! Number of rows at level 1
    integer(i4)                                          :: ncols1   ! Number of cols at level 1
    logical, dimension(:,:), allocatable                 :: mask1    ! Mask at Level 1
    real(dp)                                             :: xllcorner0, yllcorner0
    real(dp)                                             :: cellsize0
    !
    ! Dummy Variables
    integer(i4)                                          :: ii
    integer(i4), dimension(:,:),   allocatable           :: dummyI2  ! dummy, 2 dimension I4
    integer(i4), dimension(:,:),   allocatable           :: dummyI22 ! 2nd dummy, 2 dimension I4
    real(dp),    dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension DP 

    ! local variables
    character(256) :: Fname

    ! read config
    Fname = trim(InPath) // trim(num2str(iBasin, '(i3.3)')) // '_config.nc' ! '_restart.nc'
    call message('    Reading config from     ', trim(adjustl(Fname)),' ...')
 
    !
    ! level-0 information
    call get_basin_info( iBasin, 0, nrows0, ncols0, iStart= iStart0, iEnd=iEnd0, mask=mask0, &
         xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0  )
    !
    if (iBasin .eq. 1) then
       allocate( level1%nrows        (nBasins) )
       allocate( level1%ncols        (nBasins) )
       allocate( level1%xllcorner    (nBasins) )
       allocate( level1%yllcorner    (nBasins) )
       allocate( level1%cellsize     (nBasins) )
       allocate( level1%nodata_value (nBasins) )
    end if
    ! grid properties
    call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp,         &
         resolutionHydrology(iBasin) , &
         level1%nrows(iBasin), level1%ncols(iBasin), level1%xllcorner(iBasin), &
         level1%yllcorner(iBasin), level1%cellsize(iBasin), level1%nodata_value(iBasin) )
    !
    ! level-1 information
    call get_basin_info( iBasin, 1, nrows1, ncols1 )

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L0 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! check whether L0 data should be read
    if ( iBasin .eq. 1 ) then
       ! read L0 cellCoor
       allocate( dummyI2( nrows0, ncols0 ) )
       allocate( dummyI22( count(mask0), 2 ))
       call Get_NcVar( Fname, 'L0_rowCoor', dummyI2 )
       dummyI22(:,1) = pack( dummyI2, mask0 )
       call Get_NcVar( Fname, 'L0_colCoor', dummyI2 )
       dummyI22(:,2) = pack( dummyI2, mask0 )
       call append( L0_cellCoor, dummyI22)
       deallocate( dummyI22 )
       !
       call Get_NcVar( Fname, 'L0_Id', dummyI2)
       call append( L0_Id, pack(dummyI2, mask0) )
       deallocate( dummyI2 )
       !
       allocate( dummyD2( nrows0, ncols0 ) )
       call Get_NcVar( Fname, 'L0_areaCell', dummyD2 )
       call append( L0_areaCell, pack(dummyD2, mask0) )
       !
       call Get_NcVar( Fname, 'L0_slope_emp', dummyD2 )
       call append( L0_slope_emp, pack(dummyD2, mask0) )
       deallocate( dummyD2 )
    else
       if ( L0_Basin(iBasin) .ne. L0_Basin(iBasin - 1) ) then
          ! read L0 cellCoor
          allocate( dummyI2( nrows0, ncols0 ) )
          allocate( dummyI22( count(mask0), 2 ))
          call Get_NcVar( Fname, 'L0_rowCoor', dummyI2 )
          dummyI22(:,1) = pack( dummyI2, mask0 )
          call Get_NcVar( Fname, 'L0_colCoor', dummyI2 )
          dummyI22(:,2) = pack( dummyI2, mask0 )
          call append( L0_cellCoor, dummyI22)
          deallocate( dummyI22 )
          !
          call Get_NcVar( Fname, 'L0_Id', dummyI2)
          call append( L0_Id, pack(dummyI2, mask0) )
          deallocate( dummyI2 )
          !
          allocate( dummyD2( nrows0, ncols0 ) )
          call Get_NcVar( Fname, 'L0_areaCell', dummyD2 )
          call append( L0_areaCell, pack(dummyD2, mask0) )
          !
          call Get_NcVar( Fname, 'L0_slope_emp', dummyD2 )
          call append( L0_slope_emp, pack(dummyD2, mask0) )
          deallocate( dummyD2 )
       end if
    end if

    ! update L0_nCells
    L0_nCells = size(L0_Id,1)

    !------------------------------------------------------
    ! Assign whether a given soil type is present or not
    !------------------------------------------------------
    if ( iBasin .eq. 1 ) then
       allocate( soilId_isPresent(nSoilTypes) )
       soilId_isPresent(:) = 0
    end if
    !------------------------------------------------------
    ! set soil types when mpr should be performed
    !------------------------------------------------------
    if ( perform_mpr ) then
       do ii = iStart0, iEnd0
          soilId_isPresent(L0_soilId(ii)) = 1
       end do
    end if
    !
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L1 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L1 mask
    allocate( dummyI2( nrows1, ncols1 ), mask1( nrows1, ncols1) )
    call Get_NcVar( Fname,  'L1_basin_Mask', dummyI2 )
    mask1 = (dummyI2 .eq. 1_i4)
    call append( basin%L1_Mask, reshape( mask1, (/nrows1*ncols1/)))
    !
    ! update basin database
    if( iBasin .eq. 1 ) then

       allocate(basin%L1_iStart(nBasins))
       allocate(basin%L1_iEnd  (nBasins))
       allocate(basin%L1_iStartMask(nBasins))
       allocate(basin%L1_iEndMask   (nBasins))    

       ! basin information
       basin%L1_iStart(iBasin) = 1
       basin%L1_iEnd  (iBasin) = basin%L1_iStart(iBasin) + count(mask1) - 1

       basin%L1_iStartMask(iBasin) = 1
       basin%L1_iEndMask  (iBasin) = basin%L1_iStartMask(iBasin) + nrows1*ncols1 - 1

    else

       ! basin information
       basin%L1_iStart(iBasin) = basin%L1_iEnd(iBasin-1) + 1
       basin%L1_iEnd  (iBasin) = basin%L1_iStart(iBasin) + count(mask1) - 1

       basin%L1_iStartMask(iBasin) = basin%L1_iEndMask(iBasin-1) + 1
       basin%L1_iEndMask  (iBasin) = basin%L1_iStartMask(iBasin) + nrows1*ncols1 - 1

    end if
    !
    ! L1 cell Ids
    call Get_NcVar( Fname, 'L1_Id', dummyI2)
    call append( L1_Id, pack(dummyI2, mask1) )
    ! L1 cell coordinates
    allocate( dummyI22( count(mask1), 2 ))
    call Get_NcVar( Fname, 'L1_rowCoor', dummyI2 )
    dummyI22(:,1) = pack( dummyI2, mask1 )
    call Get_NcVar( Fname, 'L1_colCoor', dummyI2 )
    dummyI22(:,2) = pack( dummyI2, mask1 )
    call append( L1_cellCoor, dummyI22)
    deallocate( dummyI22 )
    ! 
    call Get_NcVar( Fname, 'L1_upBound_L0', dummyI2)
    call append( L1_upBound_L0   , pack( dummyI2, mask1) )
    !
    call Get_NcVar( Fname, 'L1_downBound_L0', dummyI2 )
    call append( L1_downBound_L0 , pack( dummyI2, mask1)  )
    !
    call Get_NcVar( Fname, 'L1_leftBound_L0', dummyI2 )
    call append( L1_leftBound_L0 , pack( dummyI2, mask1)  )
    !
    call Get_NcVar( Fname, 'L1_rightBound_L0', dummyI2 )
    call append( L1_rightBound_L0, pack( dummyI2, mask1) )
    !
    call Get_NcVar( Fname, 'L1_nTCells_L0', dummyI2 )
    call append( L1_nTCells_L0   , pack( dummyI2, mask1)    )
    deallocate( dummyI2)
    !
    allocate( dummyD2( nrows1, ncols1 ) )
    call Get_NcVar( Fname, 'L1_areaCell', dummyD2 )
    call append( L1_areaCell     , pack( dummyD2, mask1)   )

    L1_nCells = size( L1_Id, 1 )

  end subroutine read_restart_config

  ! ------------------------------------------------------------------

  !      NAME
  !         read_restart_states

  !     PURPOSE
  !>        \brief reads fluxes and state variables from file

  !>        \details read fluxes and state variables from given 
  !>        restart directory and initialises all state variables
  !>        that are initialized in the subroutine initialise,
  !>        contained in module mo_startup.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"        number of basin
  !>        \param[in] "character(256) :: InPath"        Input Path including trailing slash

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN

  !     RESTRICTIONS 
  !>        \note Restart Files must have the format, as if
  !>        it would have been written by subroutine write_restart_files

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Apr 2013
  !         Modified   R. Kumar, J. Mai    Sep. 2013  - Splitting allocation and initialization of arrays
  !         Modified   Matthias Zink       Nov. 2014  - added PET related parameter read in

  subroutine read_restart_states( iBasin, InPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_ncread,           only: Get_NcVar
    use mo_mhm_constants,    only: nRoutingStates, YearMonths_i4
    use mo_global_variables, only: processMatrix, &
         L1_fSealed, &
         L1_fForest, &
         L1_fPerm, &
         L1_Inter, &
         L1_snowPack, &
         L1_sealSTW, &
         L1_soilMoist, &
         L1_unsatSTW, &
         L1_satSTW, &
         L1_aETSoil, &
         L1_aETCanopy, &
         L1_aETSealed, &
         L1_baseflow, &
         L1_infilSoil, &
         L1_fastRunoff, &
         L1_melt, &
         L1_percol, &
         L1_preEffect, &
         L1_rain, &
         L1_runoffSeal, &
         L1_slowRunoff, &
         L1_snow, &
         L1_Throughfall, &
         L1_total_runoff, &
         L1_alpha, &
         L1_degDayInc, &
         L1_degDayMax, &
         L1_degDayNoPre, &
         L1_degDay, &
         L1_karstLoss, &
         L1_fAsp, &
         L1_HarSamCoeff, &
         L1_PrieTayAlpha, &
         L1_aeroResist, &
         L1_surfResist, &
         L1_fRoots, &
         L1_maxInter, &
         L1_kfastFlow, &
         L1_kSlowFlow, &
         L1_kBaseFlow, &
         L1_kPerco, &
         L1_soilMoistFC, &
         L1_soilMoistSat, &
         L1_soilMoistExp, &
         L1_tempThresh, &
         L1_unsatThresh, &
         L1_sealedThresh, &
         L1_wiltingPoint, &
         L11_Qmod, &
         L11_qOUT, &
         L11_qTIN, &
         L11_qTR, &
         L11_K, &
         L11_xi, &
         L11_C1, &
         L11_C2, &
         L11_FracFPimp, &
         nSoilHorizons_mHM

    implicit none

    integer(i4),    intent(in) :: iBasin
    character(256), intent(in) :: InPath ! list of Output paths per Basin

    character(256)                                    :: Fname
    integer(i4)                                       :: ii       ! loop index
    integer(i4)                                       :: s1       ! start index at level 1
    integer(i4)                                       :: e1       ! end index at level 1
    integer(i4)                                       :: ncols1   ! number of colums at level 1
    integer(i4)                                       :: nrows1   ! number of rows at level 1
    integer(i4)                                       :: ncells1  ! number of cells at level 1
    logical, dimension(:,:), allocatable              :: mask1    ! mask at level 1
    integer(i4)                                       :: s11      ! start index at level 11
    integer(i4)                                       :: e11      ! end index at level 11
    integer(i4)                                       :: ncols11  ! number of colums at level 11
    integer(i4)                                       :: nrows11  ! number of rows at level 11
    integer(i4)                                       :: ncells11 ! number of cells at level 11
    logical, dimension(:,:), allocatable              :: mask11   ! mask at level 11

    real(dp), dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension
    real(dp), dimension(:,:,:), allocatable           :: dummyD3  ! dummy, 3 dimension

    Fname = trim(InPath) // trim(num2str(iBasin, '(i3.3)')) // '_states.nc'! '_restart.nc'
    call message('    Reading states from ', trim(adjustl(Fname)),' ...')

    ! get basin information at level 1
    call get_basin_info( iBasin, 1, nrows1, ncols1, ncells=ncells1, &
         iStart=s1, iEnd=e1, mask=mask1 )

    !-------------------------------------------
    ! LAND COVER variables
    !-------------------------------------------
    allocate( dummyD2( nrows1, ncols1 ) )

    call Get_NcVar( Fname,  'L1_fSealed', dummyD2 )
    L1_fSealed(s1:e1) = pack( dummyD2, mask1 )

    call Get_NcVar( Fname,  'L1_fForest', dummyD2 )
    L1_fForest(s1:e1) = pack( dummyD2, mask1 )

    call Get_NcVar( Fname,  'L1_fPerm', dummyD2 )
    L1_fPerm(s1:e1) = pack( dummyD2, mask1 )

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------

    ! Interception
    call Get_NcVar( Fname,  'L1_Inter', dummyD2 )
    L1_inter(s1:e1) = pack( dummyD2, mask1 )

    ! Snowpack
    call Get_NcVar( Fname,  'L1_snowPack', dummyD2 )
    L1_snowPack(s1:e1) = pack( dummyD2, mask1 )

    ! Retention storage of impervious areas
    call Get_NcVar( Fname,  'L1_sealSTW', dummyD2 )
    L1_sealSTW(s1:e1) = pack( dummyD2, mask1 )

    ! upper soil storage
    call Get_NcVar( Fname,  'L1_unsatSTW', dummyD2 )
    L1_unsatSTW(s1:e1) = pack( dummyD2, mask1 )

    ! groundwater storage
    call Get_NcVar( Fname,  'L1_satSTW', dummyD2 )
    L1_satSTW(s1:e1) = pack( dummyD2, mask1 )

    ! Soil moisture of each horizon
    deallocate( dummyD2 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1, nSoilHorizons_mHM ) )

    call Get_NcVar( Fname,  'L1_soilMoist', dummyD3 )

    do ii = 1, nSoilHorizons_mHM
       L1_soilMoist(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------   

    !  soil actual ET
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_aETSoil', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_aETSoil(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    deallocate( dummyD2 )
    allocate( dummyD2( nrows1, ncols1 ) )

    ! canopy actual ET
    call Get_NcVar( Fname,  'L1_aETCanopy', dummyD2 )
    L1_aETCanopy(s1:e1) = pack( dummyD2, mask1 ) 

    ! sealed area actual ET
    call Get_NcVar( Fname,  'L1_aETSealed', dummyD2 )
    L1_aETSealed(s1:e1) = pack( dummyD2, mask1 ) 

    ! baseflow
    call Get_NcVar( Fname,  'L1_baseflow', dummyD2 )
    L1_baseflow(s1:e1) = pack( dummyD2, mask1 ) 

    ! soil in-exfiltration
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_infilSoil', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_infilSoil(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    deallocate( dummyD2 )
    allocate( dummyD2( nrows1, ncols1 ) )

    ! fast runoff
    call Get_NcVar( Fname,  'L1_fastRunoff', dummyD2 )
    L1_fastRunoff(s1:e1) = pack( dummyD2, mask1 )

    ! snow melt
    call Get_NcVar( Fname,  'L1_melt', dummyD2 )
    L1_melt(s1:e1) = pack( dummyD2, mask1 )

    ! percolation
    call Get_NcVar( Fname,  'L1_percol', dummyD2 )
    L1_percol(s1:e1) = pack( dummyD2, mask1 ) 

    ! effective precip. depth (snow melt + rain)
    call Get_NcVar( Fname,  'L1_preEffect', dummyD2 )
    L1_preEffect(s1:e1) = pack( dummyD2, mask1 )

    ! rain (liquid water)
    call Get_NcVar( Fname,  'L1_rain', dummyD2 )
    L1_rain(s1:e1) = pack( dummyD2, mask1 ) 

    ! runoff from impervious area
    call Get_NcVar( Fname,  'L1_runoffSeal', dummyD2 )
    L1_runoffSeal(s1:e1) = pack( dummyD2, mask1 )

    ! slow runoff
    call Get_NcVar( Fname,  'L1_slowRunoff', dummyD2 )
    L1_slowRunoff(s1:e1) = pack( dummyD2, mask1 ) 

    ! snow (solid water)
    call Get_NcVar( Fname,  'L1_snow', dummyD2 )
    L1_snow(s1:e1) = pack( dummyD2, mask1 )

    ! throughfall 
    call Get_NcVar( Fname,  'L1_Throughfall', dummyD2 )
    L1_Throughfall(s1:e1) = pack( dummyD2, mask1 )

    ! total runoff
    call Get_NcVar( Fname,  'L1_total_runoff', dummyD2 )
    L1_total_runoff(s1:e1) = pack( dummyD2, mask1 )

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------

    ! exponent for the upper reservoir
    call Get_NcVar( Fname,  'L1_alpha', dummyD2 )
    L1_alpha(s1:e1) = pack( dummyD2, mask1 ) 

    ! increase of the Degree-day factor per mm of increase in precipitation
    call Get_NcVar( Fname,  'L1_degDayInc', dummyD2 )
    L1_degDayInc(s1:e1) = pack( dummyD2, mask1 ) 

    ! maximum degree-day factor 
    call Get_NcVar( Fname,  'L1_degDayMax', dummyD2 )
    L1_degDayMax(s1:e1) = pack( dummyD2, mask1 ) 

    ! degree-day factor with no precipitation
    call Get_NcVar( Fname,  'L1_degDayNoPre', dummyD2 )
    L1_degDayNoPre(s1:e1) = pack( dummyD2, mask1 ) 

    ! degree-day factor
    call Get_NcVar( Fname,  'L1_degDay', dummyD2 )
    L1_degDay(s1:e1) = pack( dummyD2, mask1 ) 

    ! Karstic percolation loss
    call Get_NcVar( Fname,  'L1_karstLoss', dummyD2 )
    L1_karstLoss(s1:e1) = pack( dummyD2, mask1 ) 

    ! Fraction of roots in soil horizons    
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_fRoots', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_fRoots(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Maximum interception 
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD2( nrows1, ncols1 ) )

    call Get_NcVar( Fname, 'L1_maxInter', dummyD2 )
    L1_maxInter(s1:e1) = pack(dummyD2, mask1) 

    ! fast interflow recession coefficient 
    call Get_NcVar( Fname,  'L1_kfastFlow', dummyD2 )
    L1_kfastFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! slow interflow recession coefficient 
    call Get_NcVar( Fname,  'L1_kSlowFlow', dummyD2 )
    L1_kSlowFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! baseflow recession coefficient 
    call Get_NcVar( Fname,  'L1_kBaseFlow', dummyD2 )
    L1_kBaseFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! percolation coefficient
    call Get_NcVar( Fname,  'L1_kPerco', dummyD2 )
    L1_kPerco(s1:e1) = pack( dummyD2, mask1 ) 

    ! Soil moisture below which actual ET is reduced linearly till PWP
    deallocate( dummyD2 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_soilMoistFC', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistFC(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Saturation soil moisture for each horizon [mm]
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_soilMoistSat', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistSat(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Exponential parameter to how non-linear is the soil water retention
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_soilMoistExp', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistExp(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    deallocate( dummyD2 )
    allocate( dummyD2( nrows1, ncols1 ) )

    ! Threshold temperature for snow/rain 
    call Get_NcVar( Fname,  'L1_tempThresh', dummyD2 )
    L1_tempThresh(s1:e1) = pack( dummyD2, mask1 ) 

    ! Threshhold water depth controlling fast interflow
    call Get_NcVar( Fname,  'L1_unsatThresh', dummyD2 )
    L1_unsatThresh(s1:e1) = pack( dummyD2, mask1 )

    ! Threshhold water depth for surface runoff in sealed surfaces
    call Get_NcVar( Fname,  'L1_sealedThresh', dummyD2 )
    L1_sealedThresh(s1:e1) = pack( dummyD2, mask1 ) 

    ! Permanent wilting point
    deallocate( dummyD2, dummyD3 )
    allocate( dummyD3( nrows1, ncols1, nSoilHorizons_mHM ) )
    allocate( dummyD2( ncells1,  nSoilHorizons_mHM ) )

    call Get_NcVar( Fname, 'L1_wiltingPoint', dummyD3 )
    do ii = 1, nSoilHorizons_mHM
       L1_wiltingPoint(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    deallocate( dummyD2, dummyD3 )

    ! different parameters dependent on PET formulation
    select case (processMatrix(5,1))
    case(0) ! PET is input
       allocate( dummyD2( nrows1, ncols1 ) )

       ! PET correction factor due to terrain aspect
       call Get_NcVar( Fname,  'L1_fAsp', dummyD2 )
       L1_fAsp(s1:e1) = pack( dummyD2, mask1 ) 

       deallocate( dummyD2)
    case(1) ! HarSam
       allocate( dummyD2( nrows1, ncols1 ) )

       ! PET correction factor due to terrain aspect
       call Get_NcVar( Fname,  'L1_fAsp', dummyD2 )
       L1_fAsp(s1:e1) = pack( dummyD2, mask1 ) 

       ! Hargreaves Samani coeffiecient
       call Get_NcVar( Fname,  'L1_HarSamCoeff', dummyD2 )
       L1_HarSamCoeff(s1:e1) = pack( dummyD2, mask1 ) 

       deallocate( dummyD2)
    case(2) ! PrieTay
       allocate( dummyD3( nrows1, ncols1, YearMonths_i4) )

       ! Priestley Taylor coeffiecient (alpha)
       call Get_NcVar( Fname,  'L1_PrieTayAlpha', dummyD3 )
       do ii = 1, YearMonths_i4
          L1_PrieTayAlpha(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
       end do

       deallocate( dummyD3)
    case(3) ! PenMon
       allocate( dummyD3( nrows1, ncols1, YearMonths_i4) )

       ! aerodynamical resitance
       call Get_NcVar( Fname,  'L1_aeroResist', dummyD3 )
       do ii = 1, YearMonths_i4
          L1_aeroResist(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
       end do

       ! bulk surface resitance
       call Get_NcVar( Fname,  'L1_surfResist', dummyD3 )
       do ii = 1, YearMonths_i4
          L1_surfResist(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
       end do

       deallocate( dummyD3)
    end select

    !-------------------------------------------
    ! L11 ROUTING STATE VARIABLES, FLUXES AND
    !             PARAMETERS
    !-------------------------------------------
    if (  processMatrix(8, 1) /= 0 ) then

       ! level-11 information
       call get_basin_info( iBasin, 11, nrows11, ncols11, ncells=ncells11, &
            iStart=s11, iEnd=e11, mask=mask11 )
       allocate( dummyD2( nrows11, ncols11 ) )

       ! simulated discharge at each node
       call Get_NcVar( Fname,  'L11_Qmod', dummyD2 )
       L11_Qmod(s11:e11) = pack( dummyD2, mask11 ) 

       ! Total outflow from cells L11 at time tt
       call Get_NcVar( Fname,  'L11_qOUT', dummyD2 )
       L11_qOUT(s11:e11) = pack( dummyD2, mask11 )

       ! Total discharge inputs at t-1 and t
       deallocate( dummyD2 )
       allocate( dummyD3( nrows11, ncols11, nRoutingStates ) )
       call Get_NcVar( Fname, 'L11_qTIN', dummyD3 )
       do ii = 1, nRoutingStates
          L11_qTIN(s11:e11,ii) = pack( dummyD3(:,:,ii), mask11 )
       end do

       !  Routed outflow leaving a node
       deallocate( dummyD3 )
       allocate( dummyD3( nrows11, ncols11, nRoutingStates ) )

       allocate( dummyD2( ncells11,  nRoutingStates ) )

       call Get_NcVar( Fname, 'L11_qTR', dummyD3 )
       do ii = 1, nRoutingStates
          L11_qTR(s11:e11,ii) = pack( dummyD3(:,:,ii), mask11 )
       end do

       deallocate(dummyD2)
       allocate( dummyD2( nrows11, ncols11 ) )

       ! kappa: Muskingum travel time parameter.
       call Get_NcVar( Fname,  'L11_K', dummyD2 )
       L11_K(s11:e11) = pack( dummyD2, mask11 )

       !  xi:    Muskingum diffusion parameter
       call Get_NcVar( Fname,  'L11_xi', dummyD2 )
       L11_xi(s11:e11) = pack( dummyD2, mask11 )

       ! Routing parameter C1=f(K,xi, DT) (Chow, 25-41)
       call Get_NcVar( Fname,  'L11_C1', dummyD2 )
       L11_C1(s11:e11) = pack( dummyD2, mask11 )

       ! Routing parameter C2 =f(K,xi, DT) (Chow, 25-41)
       call Get_NcVar( Fname,  'L11_C2', dummyD2 )
       L11_C2(s11:e11) = pack( dummyD2, mask11 )

       ! Fraction of the flood plain with impervious cover
       call Get_NcVar( Fname,  'L11_FracFPimp', dummyD2 )
       L11_FracFPimp(s11:e11) = pack( dummyD2, mask11 )

       ! free memory
       deallocate( dummyD2, dummyD3 )

    end if

  end subroutine read_restart_states

END MODULE mo_restart
