!> \file mo_mpr_eval.f90
!> \brief \copybrief mo_mpr_eval
!> \details \copydetails mo_mpr_eval

!> \brief Runs MPR
!> \details Runs MPR and writes to global effective parameters
!> \authors Robert Schweppe
!> \date Feb 2018
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mpr
MODULE mo_mpr_eval

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mpr_eval

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        mpr_eval

  !    PURPOSE
  !>       \brief Runs MPR and writes to global effective parameters

  !>       \details Runs MPR and writes to global effective parameters

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), dimension(:), optional :: parameterset" a set of global parameter (gamma) to run mHM,
  !>       DIMENSION [no. of global_Parameters]

  !    HISTORY
  !>       \authors Juliane Mai, Rohini Kumar

  !>       \date Feb 2013

  ! Modifications:
  ! R. Kumar             Jun 2013 - restart_flag_states_read is passed to mhm call for the soil moisture initalisation
  ! R. Kumar             Jun 2013 - frac_sealed_city_area is added
  ! R. Kumar & S. Thober Aug 2013 - code change to incorporate output timestep during writing of the netcdf file
  ! R. Kumar             Aug 2013 - added iFlag_LAI_data_format to handle LAI options, and changed within the code made accordingly
  ! R. Kumar, J. Mai     Sep 2013 - Splitting allocation and initialization of arrays
  ! R. Kumar             Nov 2013 - update intent variables in documentation
  ! L. Samaniego         Nov 2013 - relational statements == to .eq., etc.
  ! M. Zink              Feb 2014 - added PET calculation: Hargreaves-Samani (Process 5)
  ! M. Zink              Mar 2014 - added inflow from upstream areas
  ! Stephan Thober       Jun 2014 - added chunk read for meteorological input
  ! Stephan Thober       Jun 2014 - updated flag for read_restart
  ! M. Cuntz & J. Mai    Nov 2014 - LAI input from daily, monthly or yearly files
  ! Matthias Zink        Dec 2014 - adopted inflow gauges to ignore headwater cells
  ! Stephan Thober       Aug 2015 - moved writing of daily discharge to mo_write_routing, included routing related variables from mRM
  ! David Schaefer       Aug 2015 - changed to new netcdf-writing scheme
  ! Stephan Thober       Sep 2015 - updated mrm_routing call
  ! O. Rakovec, R. Kumar Oct 2015 - added optional output for domain averaged TWS
  ! Rohini Kumar         Mar 2016 - changes for handling multiple soil database options
  ! Stephan Thober       Nov 2016 - added two options for routing
  ! Rohini Kuamr         Dec  2016 - option to handle monthly mean gridded fields of LAI
  ! Stephan Thober       Jan 2017 - added prescribed weights for tavg and pet
  ! Zink M. Demirel C.   Mar 2017 - Added Jarvis soil water stress function at SM process(3)
  ! Robert Schweppe      Dec 2017 - extracted call to mpr from inside mhm

  SUBROUTINE mpr_eval(parameterset, opti_domain_indices)

    use mo_common_variables, only : L0_LCover, l0_l1_remap, level0, level1, domainMeta
    use mo_message, only : message
    use mo_mpr_global_variables, only : L0_asp, L0_geoUnit, L0_gridded_LAI, &
                                        L0_slope_emp, L0_soilId, L1_HarSamCoeff, L1_PrieTayAlpha, L1_aeroResist, &
                                        L1_alpha, L1_degDayInc, L1_degDayMax, L1_degDayNoPre, L1_fAsp, L1_fRoots, &
                                        L1_fSealed, L1_jarvis_thresh_c1, L1_kBaseFlow, L1_kPerco, L1_kSlowFlow, &
                                        L1_karstLoss, L1_kfastFlow, L1_maxInter, L1_petLAIcorFactor, L1_sealedThresh, &
                                        L1_soilMoistExp, L1_soilMoistFC, L1_soilMoistSat, L1_surfResist, L1_tempThresh, &
                                        L1_unsatThresh, L1_wiltingPoint, &
                                        ! neutron count
                                        L1_No_Count, L1_bulkDens, L1_latticeWater, L1_COSMICL3

    use mo_multi_param_reg, only : mpr
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_clear, timer_get, timer_start, timer_stop

    implicit none

    ! a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]
    real(dp), dimension(:), intent(in), optional :: parameterset

    integer(i4), dimension(:), optional, intent(in) :: opti_domain_indices

    ! Counters
    integer(i4) :: iDomain, ii, nDomains, itimer

    ! start and end index at level 0 for current domain
    integer(i4) :: s0, e0

    ! start and end index at level 1 for current domain
    integer(i4) :: s1, e1


    !-------------------------------------------------------------------
    ! NOW call MPR
    !-------------------------------------------------------------------
    call message('  Executing MPR ...')
    itimer = 10
    call timer_start(itimer)

    !----------------------------------------
    ! loop over domains
    !----------------------------------------
    if (present(opti_domain_indices)) then
      nDomains = size(opti_domain_indices)
    else
      nDomains = domainMeta%nDomains
    end if
    do ii = 1, nDomains
      if (present(opti_domain_indices)) then
        iDomain = opti_domain_indices(ii)
      else
        iDomain = ii
      end if

      ! get domain information
      s0 = level0(domainMeta%L0DataFrom(iDomain))%iStart
      e0 = level0(domainMeta%L0DataFrom(iDomain))%iEnd
      s1 = level1(iDomain)%iStart
      e1 = level1(iDomain)%iEnd

      call mpr(level0(domainMeta%L0DataFrom(iDomain))%mask, L0_geoUnit(s0 : e0), &
              L0_soilId(s0 : e0, :), L0_asp(s0 : e0), L0_gridded_LAI(s0 : e0, :), &
              L0_LCover(s0 : e0, :), L0_slope_emp(s0 : e0), &
              pack(level0(domainMeta%L0DataFrom(iDomain))%y, level0(domainMeta%L0DataFrom(iDomain))%mask), &
              level0(domainMeta%L0DataFrom(iDomain))%Id, &
              l0_l1_remap(iDomain)%upper_bound, l0_l1_remap(iDomain)%lower_bound, &
              l0_l1_remap(iDomain)%left_bound, l0_l1_remap(iDomain)%right_bound, &
              l0_l1_remap(iDomain)%n_subcells, &
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
              L1_No_Count(s1:e1,:,:), L1_bulkDens(s1:e1,:,:), L1_latticeWater(s1:e1,:,:), L1_COSMICL3(s1:e1,:,:), &
              parameterset )

    end do
    call timer_stop(itimer)
    call message('  in ', trim(num2str(timer_get(itimer), '(F9.3)')), ' seconds.')
    call timer_clear(itimer)

  END SUBROUTINE mpr_eval

END MODULE mo_mpr_eval
