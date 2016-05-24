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
  !         None

  !     HISTORY
  !>        \author   Stephan Thober
  !>        \date     Jun 2014
  !         Modified  Matthias Zink   Nov. 2014  - added PET related parameter writing
  !                   Stephan Thober  Aug  2015  - moved write of routing states to mRM
  !                   David Schaefer  Nov  2015  - mo_netcdf
  ! ------------------------------------------------------------------ 
  subroutine write_restart_files( OutPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_init_states,      only: get_basin_info
    use mo_string_utils,     only: num2str
    use mo_netcdf,           only: NcDataset, NcDimension, NcVariable
    use mo_mhm_constants,    only: nodata_dp, nodata_i4
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
         basin, &
         L0_cellCoor    ,          & 
         L0_Id         ,           & ! Ids of grid at level-0 
         L0_slope_emp  ,           & ! Empirical quantiles of slope
         L0_areaCell,              & ! Ids of grid at level-0
         L1_areaCell,              & ! [km2] Effective area of cell at this level
         L1_Id         ,           & ! Ids of grid at level-1
         L1_cellCoor    ,          &
         L1_upBound_L0 ,           & ! Row start at finer level-0 scale 
         L1_downBound_L0,          & ! Row end at finer level-0 scale
         L1_leftBound_L0,          & ! Col start at finer level-0 scale
         L1_rightBound_L0,         & ! Col end at finer level-0 scale
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
    real(dp), dimension(:,:,:), allocatable  :: dummy_d3 ! dummy variable

    type(NcDataset)                          :: nc
    type(NcDimension)                        :: rows0, cols0, rows1, cols1, soil1, months
    type(NcVariable)                         :: var
    
    basin_loop: do iBasin = 1, size(OutPath)

       ! get Level0 information about the basin
       call get_basin_info( iBasin, 0, nrows0, ncols0, iStart=s0, iEnd=e0, mask=mask0 )

       ! get Level1 information about the basin
       call get_basin_info( iBasin, 1, nrows1, ncols1, iStart=s1, iEnd=e1, mask=mask1 )

       ! write restart file for iBasin
       Fname = trim(OutPath(iBasin)) // "mHM_restart_" // trim(num2str(iBasin, "(i3.3)")) // ".nc"
       ! print a message
       call message("    Writing Restart-file: ", trim(adjustl(Fname))," ...")

       nc     = NcDataset(fname, "w")
       rows0  = nc%setDimension("nrows0", nrows0)
       cols0  = nc%setDimension("ncols0", ncols0)
       rows1  = nc%setDimension("nrows1", nrows1)
       cols1  = nc%setDimension("ncols1", ncols1)
       soil1  = nc%setDimension("L1_soilhorizons", size( L1_soilMoist, 2))
       months = nc%setDimension("MonthsPerYear", size( L1_PrieTayAlpha, 2))

       var = nc%setVariable("L1_fSealed","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_fSealed(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","fraction of Sealed area at level 1")

       var = nc%setVariable("L1_fForest","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_fForest(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","fraction of Forest area at level 1")
      
       var = nc%setVariable("L1_fPerm","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_fPerm(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","fraction of permeable area at level 1")

       var = nc%setVariable("L1_Inter","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_inter(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Interception storage at level 1")

       var = nc%setVariable("L1_snowPack","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_snowPack(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Snowpack at level 1")

       var = nc%setVariable("L1_sealSTW","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_sealSTW(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Retention storage of impervious areas at level 1")

       allocate( dummy_d3( nrows1, ncols1, size( L1_soilMoist, 2) ) )
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_soilMoist(s1:e1,ii), mask1, nodata_dp )
       end do

       var = nc%setVariable("L1_soilMoist","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil moisture at level 1")

       var = nc%setVariable("L1_unsatSTW","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_unsatSTW(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","upper soil storage at level 1")

       var = nc%setVariable("L1_satSTW","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_satSTW(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","groundwater storage at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_aETSoil(s1:e1,ii), mask1, nodata_dp )
       end do

       var = nc%setVariable("L1_aETSoil","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil actual ET at level 1")

       var = nc%setVariable("L1_aETCanopy","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_aETCanopy(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","canopy actual ET at level 1")

       var = nc%setVariable("L1_aETSealed","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_aETSealed(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","sealed actual ET at level 1")

       var = nc%setVariable("L1_baseflow","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_baseflow(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","baseflow at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_infilSoil(s1:e1,ii), mask1, nodata_dp )
       end do

       var = nc%setVariable("L1_infilSoil","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil in-exfiltration at level 1")

       var = nc%setVariable("L1_fastRunoff","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_fastRunoff(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","fast runoff")

       var = nc%setVariable("L1_percol","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_percol(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","percolation at level 1")

       var = nc%setVariable("L1_melt","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_melt(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","snow melt at level 1")

       var = nc%setVariable("L1_preEffect","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_preEffect(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","effective precip. depth (snow melt + rain) at level 1")

       var = nc%setVariable("L1_rain","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_rain(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","rain (liquid water) at level 1")
       
       var = nc%setVariable("L1_runoffSeal","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_runoffSeal(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","runoff from impervious area at level 1")

       var = nc%setVariable("L1_slowRunoff","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_slowRunoff(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","slow runoff at level 1")

       var = nc%setVariable("L1_snow","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_snow(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","snow (solid water) at level 1")

       var = nc%setVariable("L1_Throughfall","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_Throughfall(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","throughfall at level 1")

       var = nc%setVariable("L1_total_runoff","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_total_runoff(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","total runoff at level 1")

       !-------------------------------------------
       ! EFFECTIVE PARAMETERS
       !-------------------------------------------
       var = nc%setVariable("L1_alpha","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_alpha(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","exponent for the upper reservoir at level 1")

       var = nc%setVariable("L1_degDayInc","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_degDayInc(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","increase of the Degree-day factor per mm of increase in precipitation at level 1")

       var = nc%setVariable("L1_degDayMax","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_degDayMax(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","maximum degree-day factor at level 1")

       var = nc%setVariable("L1_degDayNoPre","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_degDayNoPre(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","degree-day factor with no precipitation at level 1")

       var = nc%setVariable("L1_degDay","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_degDay(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","degree-day factor at level 1")

       var = nc%setVariable("L1_karstLoss","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_karstLoss(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Karstic percolation loss at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_fRoots(s1:e1,ii), mask1, nodata_dp )
       end do
       var = nc%setVariable("L1_fRoots","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","Fraction of roots in soil horizons at level 1")

       var = nc%setVariable("L1_maxInter","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_maxInter(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Maximum interception at level 1")

       var = nc%setVariable("L1_kfastFlow","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_kfastFlow(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","fast interflow recession coefficient at level 1")

       var = nc%setVariable("L1_kSlowFlow","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_kSlowFlow(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","slow interflow recession coefficient at level 1")

       var = nc%setVariable("L1_kBaseFlow","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_kBaseFlow(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","baseflow recession coefficient at level 1")

       var = nc%setVariable("L1_kPerco","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_kPerco(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","percolation coefficient at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_soilMoistFC(s1:e1,ii), mask1, nodata_dp )
       end do
       var = nc%setVariable("L1_soilMoistFC","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","Soil moisture below which actual ET is reduced linearly till PWP at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_soilMoistSat(s1:e1,ii), mask1, nodata_dp )
       end do
       var = nc%setVariable("L1_soilMoistSat","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","Saturation soil moisture for each horizon [mm] at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_soilMoistExp(s1:e1,ii), mask1, nodata_dp )
       end do
       var = nc%setVariable("L1_soilMoistExp","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","Exponential parameter to how non-linear is the soil water retention at level 1")

       var = nc%setVariable("L1_tempThresh","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_tempThresh(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Threshold temperature for snow/rain at level 1")

       var = nc%setVariable("L1_unsatThresh","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_unsatThresh(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Threshhold water depth controlling fast interflow at level 1")

       var = nc%setVariable("L1_sealedThresh","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_sealedThresh(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Threshhold water depth for surface runoff in sealed surfaces at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_wiltingPoint(s1:e1,ii), mask1, nodata_dp )
       end do
       var = nc%setVariable("L1_wiltingPoint","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","Permanent wilting point at level 1")

       deallocate( dummy_d3 )

       select case (processMatrix(5,1))
       case(0) ! PET is input

          var = nc%setVariable("L1_fAsp","f64",(/rows1,cols1/))
          call var%setFillValue(nodata_dp)
          call var%setData(unpack(L1_fAsp(s1:e1), mask1, nodata_dp))
          call var%setAttribute("long_name","PET correction factor due to terrain aspect at level 1")

       case(1) ! Hargreaves-Samani

          var = nc%setVariable("L1_fAsp","f64",(/rows1,cols1/))
          call var%setFillValue(nodata_dp)
          call var%setData(unpack(L1_fAsp(s1:e1), mask1, nodata_dp))
          call var%setAttribute("long_name","PET correction factor due to terrain aspect at level 1")

          var = nc%setVariable("L1_HarSamCoeff","f64",(/rows1,cols1/))
          call var%setFillValue(nodata_dp)
          call var%setData(unpack(L1_HarSamCoeff(s1:e1), mask1, nodata_dp))
          call var%setAttribute("long_name","Hargreaves-Samani coefficient")

       case(2) ! Priestley-Taylor

          allocate( dummy_d3( nrows1, ncols1, size( L1_PrieTayAlpha, 2) ) )
          do ii = 1, size( dummy_d3, 3 )
             dummy_d3(:,:,ii) = unpack( L1_PrieTayAlpha(s1:e1,ii), mask1, nodata_dp )
          end do

          var = nc%setVariable("L1_PrieTayAlpha","f64",(/rows1,cols1,months/))
          call var%setFillValue(nodata_dp)
          call var%setData(dummy_d3)
          call var%setAttribute("long_name","Priestley Taylor coeffiecient (alpha)")

          deallocate( dummy_d3 )

       case(3) ! Penman-Monteith

          allocate( dummy_d3( nrows1, ncols1, size( L1_aeroResist, 2) ) )
          do ii = 1, size( dummy_d3, 3 )
             dummy_d3(:,:,ii) = unpack( L1_aeroResist(s1:e1,ii), mask1, nodata_dp )
          end do

          var = nc%setVariable("L1_aeroResist","f64",(/rows1,cols1,months/))
          call var%setFillValue(nodata_dp)
          call var%setData(dummy_d3)
          call var%setAttribute("long_name","aerodynamical resitance")

          do ii = 1, size( dummy_d3, 3 )
             dummy_d3(:,:,ii) = unpack( L1_surfResist(s1:e1,ii), mask1, nodata_dp )
          end do

          var = nc%setVariable("L1_surfResist","f64",(/rows1,cols1,months/))
          call var%setFillValue(nodata_dp)
          call var%setData(dummy_d3)
          call var%setAttribute("long_name","bulk surface resitance")

          deallocate( dummy_d3 )

       end select

       var = nc%setVariable("L0_rowCoor","i32",(/rows0,cols0/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L0_cellCoor(s0:e0,1), mask0, nodata_i4))
       call var%setAttribute("long_name","row coordinates at Level 0")

       var = nc%setVariable("L0_colCoor","i32",(/rows0,cols0/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L0_cellCoor(s0:e0,2), mask0, nodata_i4))
       call var%setAttribute("long_name","col coordinates at Level 0")

       var = nc%setVariable("L0_Id","i32",(/rows0,cols0/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L0_Id(s0:e0), mask0, nodata_i4))
       call var%setAttribute("long_name","cell IDs at level 0")
       
       var = nc%setVariable("L0_areaCell","f64",(/rows0,cols0/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L0_areaCell(s0:e0), mask0, nodata_dp))
       call var%setAttribute("long_name","Area of a cell at level-0 [m2]")

       var = nc%setVariable("L0_slope_emp","f64",(/rows0,cols0/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L0_slope_emp(s0:e0), mask0, nodata_dp))
       call var%setAttribute("long_name","Empirical quantiles of slope")

       var = nc%setVariable("L1_basin_Mask","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(merge(1_i4, 0_i4,  &
            reshape(basin%L1_Mask(basin%L1_iStartMask(iBasin):basin%L1_iEndMask(iBasin)),&
            (/nrows1,ncols1/))))
       call var%setAttribute("long_name","Mask at Level 1")

       var = nc%setVariable("L1_Id","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L1_Id(s1:e1), mask1, nodata_i4))
       call var%setAttribute("long_name","cell IDs at level 1")

       var = nc%setVariable("L1_rowCoor","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L1_cellCoor(s1:e1,1), mask1, nodata_i4))
       call var%setAttribute("long_name","row cell Coordinates at Level 1")

       var = nc%setVariable("L1_colCoor","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L1_cellCoor(s1:e1,2), mask1, nodata_i4))
       call var%setAttribute("long_name","col cell Coordinates at Level 1")

       var = nc%setVariable("L1_upBound_L0","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L1_upBound_L0(s1:e1), mask1, nodata_i4))
       call var%setAttribute("long_name","Row start at finer level-0 scale")

       var = nc%setVariable("L1_downBound_L0","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L1_downBound_L0(s1:e1), mask1, nodata_i4))
       call var%setAttribute("long_name","Row end at finer level-0 scale")

       var = nc%setVariable("L1_leftBound_L0","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L1_leftBound_L0(s1:e1), mask1, nodata_i4))
       call var%setAttribute("long_name","Col start at finer level-0 scale")

       var = nc%setVariable("L1_rightBound_L0","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L1_rightBound_L0(s1:e1), mask1, nodata_i4))
       call var%setAttribute("long_name","Col end at finer level-0 scal")

       var = nc%setVariable("L1_nTCells_L0","i32",(/rows1,cols1/))
       call var%setFillValue(nodata_i4)
       call var%setData(unpack(L1_nTCells_L0(s1:e1), mask1, nodata_i4))
       call var%setAttribute("long_name","Total number of valid L0 cells in a given L1 cell")

       var = nc%setVariable("L1_areaCell","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_areaCell(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","Effective area of cell at this level [km2]")

       call nc%close()
       
    end do basin_loop
    
  end subroutine write_restart_files

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
  !         Modified  David Schaefer   Nov 2015 - mo_netcdf
  subroutine read_restart_config( iBasin, soilId_isPresent, InPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_append,           only: append
    use mo_netcdf,           only: NcDataset, NcVariable
    use mo_mhm_constants,    only: nodata_dp
    use mo_init_states,      only: calculate_grid_properties
    use mo_global_variables, only: L0_Basin, & ! check whether L0_Basin should be read
         perform_mpr,       & ! switch that controls whether mpr is performed or not
         L0_soilId,         & ! soil IDs at lower level
         iFlag_soilDB,      & ! options to handle different types of soil databases
         nSoilHorizons_mHM, & ! soil horizons info for mHM         
         L0_cellCoor      , & 
         L0_Id            , & ! Ids of grid at level-0 
         L0_slope_emp     , & ! Empirical quantiles of slope
         basin,             & 
         nBasins,           &
         level1,            &
         L0_nCells,         &
         L0_areaCell,       & ! Ids of grid at level-0
         L1_areaCell,       & ! [km2] Effective area of cell at this level
         nSoilTypes,        &
         resolutionHydrology, &
         L1_nCells,           &
         L1_Id         ,      & ! Ids of grid at level-1
         L1_cellCoor   ,      &
         L1_upBound_L0 ,      & ! Row start at finer level-0 scale 
         L1_downBound_L0,     & ! Row end at finer level-0 scale
         L1_leftBound_L0,     & ! Col start at finer level-0 scale
         L1_rightBound_L0,    & ! Col end at finer level-0 scale
         L1_nTCells_L0          ! Total number of valid L0 cells in a given L1 cell

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
    integer(i4)                                          :: ii, kk, nH
    integer(i4), dimension(:,:),   allocatable           :: dummyI2  ! dummy, 2 dimension I4
    integer(i4), dimension(:,:),   allocatable           :: dummyI22 ! 2nd dummy, 2 dimension I4
    real(dp),    dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension DP 

    ! local variables
    character(256)   :: Fname
    type(NcDataset)  :: nc
    type(NcVariable) :: var
    
    ! read config
    Fname = trim(InPath) // 'mHM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc' ! '_restart.nc'
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

    nc = NcDataset(fname,"r")
    ! check whether L0 data should be read
    if ( iBasin .eq. 1 ) then
       ! read L0 cellCoor
       allocate( dummyI22( count(mask0), 2 ))

       var = nc%getVariable("L0_rowCoor")
       call var%getData(dummyI2)
       dummyI22(:,1) = pack( dummyI2, mask0 )

       var = nc%getVariable("L0_colCoor")
       call var%getData(dummyI2)
       dummyI22(:,2) = pack( dummyI2, mask0 )

       call append( L0_cellCoor, dummyI22)
       deallocate( dummyI22 )
       !
       var = nc%getVariable("L0_Id")
       call var%getData(dummyI2)
       call append( L0_Id, pack(dummyI2, mask0) )
       !
       var = nc%getVariable("L0_areaCell")
       call var%getData(dummyD2)
       call append( L0_areaCell, pack(dummyD2, mask0) )
       !
       var = nc%getVariable("L0_slope_emp")
       call var%getData(dummyD2)
       call append( L0_slope_emp, pack(dummyD2, mask0) )
    else
       if ( L0_Basin(iBasin) .ne. L0_Basin(iBasin - 1) ) then
          ! read L0 cellCoor
          allocate( dummyI22( count(mask0), 2 ))

          var = nc%getVariable("L0_rowCoor")
          call var%getData(dummyI2)
          dummyI22(:,1) = pack( dummyI2, mask0 )

          var = nc%getVariable("L0_colCoor")
          call var%getData(dummyI2)
          dummyI22(:,2) = pack( dummyI2, mask0 )

          call append( L0_cellCoor, dummyI22)
          deallocate( dummyI22 )
          !
          var = nc%getVariable("L0_Id")
          call var%getData(dummyI2)
          call append( L0_Id, pack(dummyI2, mask0) )
          !
          var = nc%getVariable("L0_areaCell")
          call var%getData(dummyD2)
          call append( L0_areaCell, pack(dummyD2, mask0) )
          !
          var = nc%getVariable("L0_slope_emp")
          call var%getData(dummyD2)
          call append( L0_slope_emp, pack(dummyD2, mask0) )
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
       nH = 1 !> by default; when iFlag_soilDB = 0
       if ( iFlag_soilDB .eq. 1 ) nH = nSoilHorizons_mHM
       do ii = 1, nH
          do kk = iStart0, iEnd0
             soilId_isPresent( L0_soilId(kk,ii) ) = 1
          end do
       end do
    end if
    !
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! Read L1 variables <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ! read L1 mask
    allocate( mask1( nrows1, ncols1) )
    var = nc%getVariable("L1_basin_Mask")
    call var%getData(dummyI2)
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
    var = nc%getVariable("L1_Id")
    call var%getData(dummyI2)
    call append( L1_Id, pack(dummyI2, mask1) )

    ! L1 cell coordinates
    allocate( dummyI22( count(mask1), 2 ))
    var = nc%getVariable("L1_rowCoor")
    call var%getData(dummyI2)
    dummyI22(:,1) = pack( dummyI2, mask1 )
    var = nc%getVariable("L1_colCoor")
    call var%getData(dummyI2)
    dummyI22(:,2) = pack( dummyI2, mask1 )
    call append( L1_cellCoor, dummyI22)
    deallocate( dummyI22 )
    ! 
    var = nc%getVariable("L1_upBound_L0")
    call var%getData(dummyI2)
    call append( L1_upBound_L0   , pack( dummyI2, mask1) )
    !
    var = nc%getVariable("L1_downBound_L0")
    call var%getData(dummyI2)
    call append( L1_downBound_L0 , pack( dummyI2, mask1)  )
    !
    var = nc%getVariable("L1_leftBound_L0")
    call var%getData(dummyI2)
    call append( L1_leftBound_L0 , pack( dummyI2, mask1)  )
    !
    var = nc%getVariable("L1_rightBound_L0")
    call var%getData(dummyI2)
    call append( L1_rightBound_L0, pack( dummyI2, mask1) )
    !
    var = nc%getVariable("L1_nTCells_L0")
    call var%getData(dummyI2)
    call append( L1_nTCells_L0   , pack( dummyI2, mask1)    )
    !
    var = nc%getVariable("L1_areaCell")
    call var%getData(dummyD2)
    call append( L1_areaCell     , pack( dummyD2, mask1)   )

    L1_nCells = size( L1_Id, 1 )

    call nc%close()

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
  !         Modified   Stephan Thober      Aug  2015  - moved read of routing states to mRM
  !         Modified   David Schaefer      Nov  2015  - mo_netcdf
  subroutine read_restart_states( iBasin, InPath )

    use mo_kind,             only: i4, dp
    ! use mo_message,          only: message
    use mo_netcdf,           only: NcDataset, NcVariable
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_mhm_constants,    only: YearMonths_i4
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

    real(dp), dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension
    real(dp), dimension(:,:,:), allocatable           :: dummyD3  ! dummy, 3 dimension

    type(NcDataset) :: nc
    type(NcVariable) :: var
    
    Fname = trim(InPath) // 'mHM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc'
    ! call message('    Reading states from ', trim(adjustl(Fname)),' ...')

    
    ! get basin information at level 1
    call get_basin_info( iBasin, 1, nrows1, ncols1, ncells=ncells1, &
         iStart=s1, iEnd=e1, mask=mask1 )

    nc = NcDataset(fname,"r")

    !-------------------------------------------
    ! LAND COVER variables
    !-------------------------------------------

    var = nc%getVariable("L1_fSealed")
    call var%getData(dummyD2)
    L1_fSealed(s1:e1) = pack( dummyD2, mask1 )

    var = nc%getVariable("L1_fForest")
    call var%getData(dummyD2)
    L1_fForest(s1:e1) = pack( dummyD2, mask1 )

    var = nc%getVariable("L1_fPerm")
    call var%getData(dummyD2)
    L1_fPerm(s1:e1) = pack( dummyD2, mask1 )

    !-------------------------------------------
    ! STATE VARIABLES
    !-------------------------------------------

    ! Interception
    var = nc%getVariable("L1_Inter")
    call var%getData(dummyD2)
    L1_inter(s1:e1) = pack( dummyD2, mask1 )

    ! Snowpack
    var = nc%getVariable("L1_snowPack")
    call var%getData(dummyD2)
    L1_snowPack(s1:e1) = pack( dummyD2, mask1 )

    ! Retention storage of impervious areas
    var = nc%getVariable("L1_sealSTW")
    call var%getData(dummyD2)
    L1_sealSTW(s1:e1) = pack( dummyD2, mask1 )

    ! upper soil storage
    var = nc%getVariable("L1_unsatSTW")
    call var%getData(dummyD2)
    L1_unsatSTW(s1:e1) = pack( dummyD2, mask1 )

    ! groundwater storage
    var = nc%getVariable("L1_satSTW")
    call var%getData(dummyD2)
    L1_satSTW(s1:e1) = pack( dummyD2, mask1 )
    
    ! Soil moisture of each horizon
    var = nc%getVariable("L1_soilMoist")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoist(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    !-------------------------------------------
    ! FLUXES
    !-------------------------------------------   

    !  soil actual ET
    var = nc%getVariable("L1_aETSoil")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_aETSoil(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! canopy actual ET
    var = nc%getVariable("L1_aETCanopy")
    call var%getData(dummyD2)
    L1_aETCanopy(s1:e1) = pack( dummyD2, mask1 ) 

    ! sealed area actual ET
    var = nc%getVariable("L1_aETSealed")
    call var%getData(dummyD2)
    L1_aETSealed(s1:e1) = pack( dummyD2, mask1 ) 

    ! baseflow
    var = nc%getVariable("L1_baseflow")
    call var%getData(dummyD2)
    L1_baseflow(s1:e1) = pack( dummyD2, mask1 ) 

    ! soil in-exfiltration
    var = nc%getVariable("L1_infilSoil")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_infilSoil(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! fast runoff
    var = nc%getVariable("L1_fastRunoff")
    call var%getData(dummyD2)
    L1_fastRunoff(s1:e1) = pack( dummyD2, mask1 )

    ! snow melt
    var = nc%getVariable("L1_melt")
    call var%getData(dummyD2)
    L1_melt(s1:e1) = pack( dummyD2, mask1 )

    ! percolation
    var = nc%getVariable("L1_percol")
    call var%getData(dummyD2)
    L1_percol(s1:e1) = pack( dummyD2, mask1 ) 

    ! effective precip. depth (snow melt + rain)
    var = nc%getVariable("L1_preEffect")
    call var%getData(dummyD2)
    L1_preEffect(s1:e1) = pack( dummyD2, mask1 )

    ! rain (liquid water)
    var = nc%getVariable("L1_rain")
    call var%getData(dummyD2)
    L1_rain(s1:e1) = pack( dummyD2, mask1 ) 

    ! runoff from impervious area
    var = nc%getVariable("L1_runoffSeal")
    call var%getData(dummyD2)
    L1_runoffSeal(s1:e1) = pack( dummyD2, mask1 )

    ! slow runoff
    var = nc%getVariable("L1_slowRunoff")
    call var%getData(dummyD2)
    L1_slowRunoff(s1:e1) = pack( dummyD2, mask1 ) 

    ! snow (solid water)
    var = nc%getVariable("L1_snow")
    call var%getData(dummyD2)
    L1_snow(s1:e1) = pack( dummyD2, mask1 )

    ! throughfall 
    var = nc%getVariable("L1_Throughfall")
    call var%getData(dummyD2)
    L1_Throughfall(s1:e1) = pack( dummyD2, mask1 )

    ! total runoff
    var = nc%getVariable("L1_total_runoff")
    call var%getData(dummyD2)
    L1_total_runoff(s1:e1) = pack( dummyD2, mask1 )

    !-------------------------------------------
    ! EFFECTIVE PARAMETERS
    !-------------------------------------------

    ! exponent for the upper reservoir
    var = nc%getVariable("L1_alpha")
    call var%getData(dummyD2)
    L1_alpha(s1:e1) = pack( dummyD2, mask1 ) 

    ! increase of the Degree-day factor per mm of increase in precipitation
    var = nc%getVariable("L1_degDayInc")
    call var%getData(dummyD2)
    L1_degDayInc(s1:e1) = pack( dummyD2, mask1 ) 

    ! maximum degree-day factor 
    var = nc%getVariable("L1_degDayMax")
    call var%getData(dummyD2)
    L1_degDayMax(s1:e1) = pack( dummyD2, mask1 ) 

    ! degree-day factor with no precipitation
    var = nc%getVariable("L1_degDayNoPre")
    call var%getData(dummyD2)
    L1_degDayNoPre(s1:e1) = pack( dummyD2, mask1 ) 

    ! degree-day factor
    var = nc%getVariable("L1_degDay")
    call var%getData(dummyD2)
    L1_degDay(s1:e1) = pack( dummyD2, mask1 ) 

    ! Karstic percolation loss
    var = nc%getVariable("L1_karstLoss")
    call var%getData(dummyD2)
    L1_karstLoss(s1:e1) = pack( dummyD2, mask1 ) 

    ! Fraction of roots in soil horizons    
    var = nc%getVariable("L1_fRoots")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_fRoots(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Maximum interception 
    var = nc%getVariable("L1_maxInter")
    call var%getData(dummyD2)
    L1_maxInter(s1:e1) = pack(dummyD2, mask1) 

    ! fast interflow recession coefficient 
    var = nc%getVariable("L1_kfastFlow")
    call var%getData(dummyD2)
    L1_kfastFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! slow interflow recession coefficient 
    var = nc%getVariable("L1_kSlowFlow")
    call var%getData(dummyD2)
    L1_kSlowFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! baseflow recession coefficient 
    var = nc%getVariable("L1_kBaseFlow")
    call var%getData(dummyD2)
    L1_kBaseFlow(s1:e1) = pack( dummyD2, mask1 ) 

    ! percolation coefficient
    var = nc%getVariable("L1_kPerco")
    call var%getData(dummyD2)
    L1_kPerco(s1:e1) = pack( dummyD2, mask1 ) 

    ! Soil moisture below which actual ET is reduced linearly till PWP

    var = nc%getVariable("L1_soilMoistFC")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistFC(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Saturation soil moisture for each horizon [mm]
    var = nc%getVariable("L1_soilMoistSat")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistSat(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Exponential parameter to how non-linear is the soil water retention
    var = nc%getVariable("L1_soilMoistExp")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_soilMoistExp(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! Threshold temperature for snow/rain 
    var = nc%getVariable("L1_tempThresh")
    call var%getData(dummyD2)
    L1_tempThresh(s1:e1) = pack( dummyD2, mask1 ) 

    ! Threshhold water depth controlling fast interflow
    var = nc%getVariable("L1_unsatThresh")
    call var%getData(dummyD2)
    L1_unsatThresh(s1:e1) = pack( dummyD2, mask1 )

    ! Threshhold water depth for surface runoff in sealed surfaces
    var = nc%getVariable("L1_sealedThresh")
    call var%getData(dummyD2)
    L1_sealedThresh(s1:e1) = pack( dummyD2, mask1 ) 

    ! Permanent wilting point
    var = nc%getVariable("L1_wiltingPoint")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_wiltingPoint(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do

    ! different parameters dependent on PET formulation
    select case (processMatrix(5,1))
    case(0) ! PET is input

       ! PET correction factor due to terrain aspect
       var = nc%getVariable("L1_fAsp")
       call var%getData(dummyD2)
       L1_fAsp(s1:e1) = pack( dummyD2, mask1 ) 

    case(1) ! Hargreaves-Samani

       ! PET correction factor due to terrain aspect
       var = nc%getVariable("L1_fAsp")
       call var%getData(dummyD2)
       L1_fAsp(s1:e1) = pack( dummyD2, mask1 ) 

       ! Hargreaves Samani coeffiecient
       var = nc%getVariable("L1_HarSamCoeff")
       call var%getData(dummyD2)
       L1_HarSamCoeff(s1:e1) = pack( dummyD2, mask1 ) 

    case(2) ! Priestely-Taylor

       ! Priestley Taylor coeffiecient (alpha)
       var = nc%getVariable("L1_PrieTayAlpha")
       call var%getData(dummyD3)
       do ii = 1, YearMonths_i4
          L1_PrieTayAlpha(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
       end do

    case(3) ! Penman-Monteith

       ! aerodynamical resitance
       var = nc%getVariable("L1_aeroResist")
       call var%getData(dummyD3)
       do ii = 1, YearMonths_i4
          L1_aeroResist(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
       end do

       ! bulk surface resitance
       var = nc%getVariable("L1_surfResist")
       call var%getData(dummyD3)
       do ii = 1, YearMonths_i4
          L1_surfResist(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
       end do

    end select

    call nc%close()

  end subroutine read_restart_states

END MODULE mo_restart
