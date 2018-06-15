!> \file mo_mpr_eval.f90

!> \brief Runs MPR and writes to global effective parameters

!> \details  Runs MPR and writes to global effective parameters

!> \authors Robert Schweppe
!> \date Feb 2018

MODULE mo_mpr_eval

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mpr_eval

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          mpr_eval

  !>        \brief Runs MPR and writes to global effective parameters

  !>        \details Runs MPR and writes to global effective parameters

  !     INTENT(IN)
  !>       \param[in] "real(dp), dimension(:)    :: parameterset"
  !>          a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]

  !     INTENT(INOUT)
  !         None


  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !>        \param[out] "real(dp), dimension(:,:), optional  :: runoff"
  !>           returns runoff time series, DIMENSION [nTimeSteps, nGaugesTotal]
  !>        \param[out] "real(dp), dimension(:,:), optional  :: sm_opti"
  !>           returns soil moisture time series for all grid cells (of multiple basins concatenated),
  !>            DIMENSION [nCells, nTimeSteps]
  !>        \param[out] "real(dp), dimension(:,:), optional  :: basin_avg_tws"
  !>           returns basin averaged total water storage time series, DIMENSION [nTimeSteps, nBasins]
  !>        \param[out] "real(dp), dimension(:,:), optional  :: neutron_opti"
  !>           returns neuton counts time series for all grid cells (of multiple basins concatenated),
  !>            DIMENSION [nCells, nTimeSteps]
  !>        \param[out] "real(dp), dimension(:,:), optional  :: et_opti"
  !>           returns evapotranspiration time series for all grid cells (of multiple basins concatenated),
  !>           DIMENSION [nCells, nTimeSteps]

  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai, Rohini Kumar
  !>        \date Feb 2013
  !         Modified, R. Kumar,             Jun 2013 - restart_flag_states_read is passed to mhm call
  !                                                    for the soil moisture initalisation
  !                   R. Kumar,             Jun 2013 - frac_sealed_city_area is added
  !                   R. Kumar & S. Thober, Aug 2013 - code change to incorporate output timestep
  !                                                    during writing of the netcdf file
  !                   R. Kumar,             Aug 2013 - added iFlag_LAI_data_format to handle LAI options,
  !                                                    and changed within the code made accordingly
  !                   R. Kumar, J. Mai,     Sep 2013 - Splitting allocation and initialization of arrays
  !                   R. Kumar              Nov 2013 - update intent variables in documentation
  !                   L. Samaniego,         Nov 2013 - relational statements == to .eq., etc.
  !                   M. Zink,              Feb 2014 - added PET calculation: Hargreaves-Samani (Process 5)
  !                   M. Zink,              Mar 2014 - added inflow from upstream areas
  !                   Stephan Thober,       Jun 2014 - added chunk read for meteorological input
  !                   Stephan Thober,       Jun 2014 - updated flag for read_restart
  !                   M. Cuntz & J. Mai,    Nov 2014 - LAI input from daily, monthly or yearly files
  !                   Matthias Zink,        Dec 2014 - adopted inflow gauges to ignore headwater cells
  !                   Stephan Thober,       Aug 2015 - moved writing of daily discharge to mo_write_routing,
  !                                                    included routing related variables from mRM
  !                   David Schaefer,       Aug 2015 - changed to new netcdf-writing scheme
  !                   Stephan Thober,       Sep 2015 - updated mrm_routing call
  !                   O. Rakovec, R. Kumar, Oct 2015 - added optional output for basin averaged TWS
  !                   Rohini Kumar,         Mar 2016 - changes for handling multiple soil database options
  !                   Stephan Thober,       Nov 2016 - added two options for routing
  !                   Rohini Kuamr,         Dec  2016 - option to handle monthly mean gridded fields of LAI
  !                   Stephan Thober,       Jan 2017 - added prescribed weights for tavg and pet
  !                   Zink M. Demirel C.,   Mar 2017 - Added Jarvis soil water stress function at SM process(3)
  !                   Robert Schweppe,      Dec 2017 - extracted call to mpr from inside mhm


  SUBROUTINE mpr_eval(parameterset)
    use mo_string_utils, only : num2str
    use mo_message, only : message
    use mo_multi_param_reg, only : mpr
    use mo_mpr_global_variables, only : &
            L0_soilId, L0_slope_emp, &
            L0_asp, L0_geoUnit, &
            L0_gridded_LAI, &
            L1_fSealed, L1_alpha, L1_degDayInc, &
            L1_degDayMax, &
            L1_degDayNoPre, L1_fAsp, L1_petLAIcorFactor, L1_HarSamCoeff, &
            L1_PrieTayAlpha, L1_aeroResist, L1_surfResist, &
            L1_fRoots, L1_maxInter, L1_karstLoss, L1_kfastFlow, &
            L1_kSlowFlow, L1_kBaseFlow, L1_kPerco, &
            L1_soilMoistFC, L1_soilMoistSat, L1_soilMoistExp, &
            L1_jarvis_thresh_c1, &
            L1_tempThresh, L1_unsatThresh, L1_sealedThresh, &
            L1_wiltingPoint
    use mo_common_variables, only : &
            L0_LCover, &
            nBasins, &
            l0_l1_remap, &
            level0, L0_Basin, &
            level1
    use mo_timer, only : timer_get, timer_start, timer_stop, timer_clear

    real(dp), dimension(:), intent(in), optional :: parameterset
    ! counters and indexes
    integer(i4) :: iBasin, itimer   ! Counters
    integer(i4) :: s0, e0           ! start and end index at level 0 for current basin
    integer(i4) :: s1, e1           ! start and end index at level 1 for current basin

    !-------------------------------------------------------------------
    ! NOW call MPR
    !-------------------------------------------------------------------
    call message('  Executing MPR ...')
    itimer = 10
    call timer_start(itimer)

    !----------------------------------------
    ! loop over basins
    !----------------------------------------
    do iBasin = 1, nBasins

      ! get basin information
      s0 = level0(L0_Basin(iBasin))%iStart
      e0 = level0(L0_Basin(iBasin))%iEnd
      s1 = level1(iBasin)%iStart
      e1 = level1(iBasin)%iEnd

      call mpr(level0(L0_Basin(iBasin))%mask, L0_geoUnit(s0 : e0), &
              L0_soilId(s0 : e0, :), L0_asp(s0 : e0), L0_gridded_LAI(s0 : e0, :), &
              L0_LCover(s0 : e0, :), L0_slope_emp(s0 : e0), &
              pack(level0(L0_Basin(iBasin))%y, level0(L0_Basin(iBasin))%mask), &
              level0(L0_Basin(iBasin))%Id, &
              l0_l1_remap(iBasin)%upper_bound, l0_l1_remap(iBasin)%lower_bound, &
              l0_l1_remap(iBasin)%left_bound, l0_l1_remap(iBasin)%right_bound, &
              l0_l1_remap(iBasin)%n_subcells, &
              L1_fSealed(s1 : e1, :, :), &
              L1_alpha(s1 : e1, :, :), L1_degDayInc(s1 : e1, :, :), L1_degDayMax(s1 : e1, :, :), &
              L1_degDayNoPre(s1 : e1, :, :), L1_fAsp(s1 : e1, :, :), L1_HarSamCoeff(s1 : e1, :, :), &
              L1_PrieTayAlpha(s1 : e1, :, :), L1_aeroResist(s1 : e1, :, :), L1_surfResist(s1 : e1, :, :), &
              L1_fRoots(s1 : e1, :, :), L1_kFastFlow(s1 : e1, :, :), &
              L1_kSlowFlow(s1 : e1, :, :), L1_kBaseFlow(s1 : e1, :, :), L1_kPerco(s1 : e1, :, :), &
              L1_karstLoss(s1 : e1, :, :), L1_soilMoistFC(s1 : e1, :, :), L1_soilMoistSat(s1 : e1, :, :), &
              L1_soilMoistExp(s1 : e1, :, :), L1_jarvis_thresh_c1(s1 : e1, :, :), &
              L1_tempThresh(s1 : e1, :, :), L1_unsatThresh(s1 : e1, :, :), L1_sealedThresh(s1 : e1, :, :), &
              L1_wiltingPoint(s1 : e1, :, :), L1_maxInter(s1 : e1, :, :), L1_petLAIcorFactor(s1 : e1, :, :), &
              parameterset)

    end do
    call timer_stop(itimer)
    call message('  in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
    call timer_clear(itimer)

  END SUBROUTINE mpr_eval

END MODULE mo_mpr_eval