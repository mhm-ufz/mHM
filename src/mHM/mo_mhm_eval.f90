!> \file mo_mhm_eval.f90
!> \brief   \copybrief mo_mhm_eval
!> \details \copydetails mo_mhm_eval

!> \brief Runs mhm with a specific parameter set and returns required variables, e.g. runoff.
!> \details Runs mhm with a specific parameter set and returns required variables, e.g. runoff.
!> \authors Juliane Mai, Rohini Kumar
!> \date Feb 2013
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_mhm_eval

  use mo_kind, only : i4, dp
  use mo_optimization_types, only : optidata_sim
  use mo_mhm_interface_run, only : &
    mhm_interface_run_prepare, &
    mhm_interface_run_get_ndomains, &
    mhm_interface_run_prepare_domain, &
    mhm_interface_run_finished, &
    mhm_interface_run_do_time_step, &
    mhm_interface_run_write_output, &
    mhm_interface_run_update_optisim, &
    mhm_interface_run_finalize_domain, &
    mhm_interface_run_finalize

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mhm_eval

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        mhm_eval

  !    PURPOSE
  !>       \brief Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

  !>       \details Runs mhm with a specific parameter set and returns required variables, e.g. runoff.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: parameterset" a set of global parameter (gamma) to run mHM, DIMENSION
  !>       [no. of global_Parameters]

  !    INTENT(OUT), OPTIONAL
  !>       \param[out] "real(dp), dimension(:, :), optional :: runoff"        returns runoff time series, DIMENSION
  !>       [nTimeSteps, nGaugesTotal]
  !>       \param[out] "real(dp), dimension(:, :), optional :: sm_opti"       returns soil moisture time series for all
  !>       grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
  !>       time series, DIMENSION [nTimeSteps, nDomains]
  !>       \param[out] "real(dp), dimension(:, :), optional :: neutrons_opti" dim1=ncells, dim2=time
  !>       \param[out] "real(dp), dimension(:, :), optional :: et_opti"       returns evapotranspiration time series for
  !>       \param[out] "real(dp), dimension(:, :), optional :: tws_opti"      returns tws time series
  !>       all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]

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
  ! O. Rakovec, R. Kumar Oct 2015 - added optional output for Domain averaged TWS
  ! Rohini Kumar         Mar 2016 - changes for handling multiple soil database options
  ! Stephan Thober       Nov 2016 - added two options for routing
  ! Rohini Kuamr         Dec 2016 - option to handle monthly mean gridded fields of LAI
  ! Stephan Thober       Jan 2017 - added prescribed weights for tavg and pet
  ! Zink M. Demirel C.   Mar 2017 - Added Jarvis soil water stress function at SM process(3)
  ! Robert Schweppe      Dec 2017 - extracted call to mpr from inside mhm
  ! Robert Schweppe      Jun 2018 - refactoring and reformatting
  ! Stephan Thober       Jan 2022 - added nTstepForcingDay and is_hourly_forcing flag

  SUBROUTINE mhm_eval(parameterset, opti_domain_indices, runoff, smOptiSim, neutronsOptiSim, etOptiSim, twsOptiSim, BFI)
    implicit none

    !> a set of global parameter (gamma) to run mHM, DIMENSION [no. of global_Parameters]
    real(dp), dimension(:), intent(in) :: parameterset
    !> selected domains for optimization
    integer(i4), dimension(:), optional, intent(in) :: opti_domain_indices
    !> returns runoff time series, DIMENSION [nTimeSteps, nGaugesTotal]
    real(dp), dimension(:, :), allocatable, optional, intent(out) :: runoff
    !> returns soil moisture time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: smOptiSim
    !> dim1=ncells, dim2=time
    type(optidata_sim), dimension(:), optional, intent(inout) :: neutronsOptiSim
    !> returns evapotranspiration time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: etOptiSim
    !> returns tws time series for all grid cells (of multiple Domains concatenated),DIMENSION [nCells, nTimeSteps]
    type(optidata_sim), dimension(:), optional, intent(inout) :: twsOptiSim
    !> baseflow index, dim1=domainID
    real(dp), dimension(:), allocatable, optional, intent(out) :: BFI

    ! number of domains simulated in this mhm_eval run. Depends on opti_function
    integer(i4) :: nDomains, ii

    ! flag to check if time loop is finished
    logical :: time_loop_finished

    ! prepare the mhm run
    call mhm_interface_run_prepare(parameterset, opti_domain_indices, present(runoff), present(BFI))

    ! get number of domains to loop over
    call mhm_interface_run_get_ndomains(nDomains)

    ! loop over Domains
    DomainLoop: do ii = 1, nDomains

      ! prepare current domain
      call mhm_interface_run_prepare_domain(ii, etOptiSim, twsOptiSim, neutronsOptiSim, smOptiSim)

      ! run time-loop at least once
      time_loop_finished = .false.

      ! Loop over time
      TimeLoop: do while(.not. time_loop_finished)

        ! do one time-step on current domain
        call mhm_interface_run_do_time_step()

        ! write output
        call mhm_interface_run_write_output()

        ! update optisim data
        call mhm_interface_run_update_optisim(etOptiSim, twsOptiSim, neutronsOptiSim, smOptiSim)

        ! check whether to run the time-loop further
        call mhm_interface_run_finished(time_loop_finished)

      end do TimeLoop !<< TIME STEPS LOOP

      ! finalize domain
      call mhm_interface_run_finalize_domain()

    end do DomainLoop !<< Domain LOOP

    ! SET RUNOFF OUTPUT VARIABLE; reset init-flag for MPR
    call mhm_interface_run_finalize(runoff, BFI)

  end SUBROUTINE mhm_eval

END MODULE mo_mhm_eval
