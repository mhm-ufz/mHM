!>       \file mo_mrm_write.f90

!>       \brief write of discharge and restart files

!>       \details This module contains the subroutines for
!>       writing the discharge files and optionally the restart
!>       files.

!>       \authors Stephan Thober

!>       \date Aug 2015

! Modifications:

module mo_mrm_write

  use mo_kind, only : i4, dp
  use mo_mrm_write_fluxes_states, only : OutputDataset

  implicit none

  public :: mrm_write
  public :: mrm_write_output_fluxes
  public :: mrm_write_optinamelist
  public :: mrm_write_optifile

  private

  ! counters for write_output_fluxes
  integer(i4) :: day_counter ! for daily output
  integer(i4) :: month_counter ! for monthly output
  integer(i4) :: year_counter ! for yearly output
  integer(i4) :: average_counter ! for averaging output
  type(OutputDataset) :: nc ! netcdf Output Dataset

contains

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_write

  !    PURPOSE
  !>       \brief write discharge and restart files

  !>       \details First, this subroutine calls the writing or restart files that only
  !>       succeeds if it happens after the write of mHM restart files because
  !>       mHM restart files must exist. Second, simulated discharge is aggregated to the daily
  !>       scale and then written to file jointly with observed discharge

  !    HISTORY
  !>       \authors Juliane Mai, Rohini Kumar & Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_write

    use mo_common_mhm_mrm_variables, only : evalPer, mrm_coupling_mode, nTstepDay, simPer, warmingDays
    use mo_common_variables, only : dirRestartOut, domainMeta, write_restart
    use mo_mrm_global_variables, only : basin_mrm, &
                                        gauge, mRM_runoff, nGaugesTotal
    use mo_mrm_restart, only : mrm_write_restart

    implicit none

    integer(i4) :: domainID, iDomain

    integer(i4) :: iDay, iS, iE

    integer(i4) :: ii

    integer(i4) :: tt

    integer(i4) :: gg

    integer(i4) :: nTimeSteps

    real(dp), dimension(:, :), allocatable :: d_Qmod


    ! --------------------------------------------------------------------------
    ! WRITE CONFIG FILE
    ! --------------------------------------------------------------------------
    if (mrm_coupling_mode .eq. 0_i4) call write_configfile()

    ! --------------------------------------------------------------------------
    ! WRITE RESTART
    ! --------------------------------------------------------------------------
    if (write_restart) then
      do iDomain = 1, domainMeta%nDomains
        domainID = domainMeta%indices(iDomain)
        call mrm_write_restart(iDomain, domainID, dirRestartOut)
      end do
    end if

    ! --------------------------------------------------------------------------
    ! STORE DAILY DISCHARGE TIMESERIES OF EACH GAUGING STATION 
    ! FOR SIMULATIONS DURING THE EVALUATION PERIOD
    !
    !  **** AT DAILY TIME STEPS ****
    ! Note:: Observed Q are stored only for the evaluation period and not for
    !        the warming days
    ! --------------------------------------------------------------------------
    ii = maxval(evalPer(1 : domainMeta%nDomains)%julEnd - evalPer(1 : domainMeta%nDomains)%julStart + 1)
    allocate(d_Qmod(ii, nGaugesTotal))
    d_Qmod = 0.0_dp

    ! loop over domains
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      nTimeSteps = (simPer(iDomain)%julEnd - simPer(iDomain)%julStart + 1) * NTSTEPDAY
      iDay = 0
      ! loop over timesteps
      do tt = warmingDays(iDomain) * NTSTEPDAY + 1, nTimeSteps, NTSTEPDAY
        iS = tt
        iE = tt + NTSTEPDAY - 1
        iDay = iDay + 1
        ! over gauges
        do gg = 1, basin_mrm(iDomain)%nGauges
          d_Qmod(iDay, basin_mrm(iDomain)%gaugeIndexList(gg)) = &
                  sum(mRM_runoff(iS : iE, basin_mrm(iDomain)%gaugeIndexList(gg))) / real(NTSTEPDAY, dp)
        end do
        !
      end do
    end do
    ! write in an ASCII file          ! OBS[nModeling_days X nGauges_total] , SIM[nModeling_days X nGauges_total] 
    call write_daily_obs_sim_discharge(gauge%Q(:, :), d_Qmod(:, :))
    ! free space
    deallocate(d_Qmod)
  end subroutine mrm_write
  ! ------------------------------------------------------------------

  !    NAME
  !        write_configfile

  !    PURPOSE
  !>       \brief This modules writes the results of the configuration into an ASCII-file

  !>       \details TODO: add description

  !    HISTORY
  !>       \authors Christoph Schneider

  !>       \date May 2013

  ! Modifications:
  ! Juliane Mai    May 2013 - module version and documentation
  ! Stephan Thober Jun 2014 - bug fix in L11 config print out
  ! Stephan Thober Jun 2014 - updated read_restart
  ! Rohini, Luis   Jul 2015 - updated version, L1 level prints
  ! Stephan Thober Sep 2015 - updated write of stream network
  ! Stephan Thober Nov 2016 - adapted write to selected case for routing process
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  Subroutine write_configfile

    use mo_common_constants, only : nodata_dp
    use mo_common_file, only : file_config, uconfig
    use mo_common_mHM_mRM_variables, only : LCyearId, SimPer, evalPer, mrm_coupling_mode, read_restart, &
                                            resolutionRouting, timeStep, warmPer
    use mo_common_variables, only : L0_Domain, LC_year_end, LC_year_start, LCfilename, &
                                    dirConfigOut, dirLCover, dirMorpho, dirOut, dirRestartOut, global_parameters, &
                                    global_parameters_name, level0, level1, domainMeta, nLCoverScene, processMatrix, &
                                    resolutionHydrology, write_restart
    use mo_kind, only : dp, i4
    use mo_message, only : message
    use mo_mrm_file, only : version
    use mo_mrm_global_variables, only : InflowGauge, L11_L1_Id, L11_fromN, L11_label, &
                                        L11_length, L11_netPerm, L11_rOrder, L11_slope, L11_toN, L1_L11_Id, basin_mrm, &
                                        dirGauges, dirTotalRunoff, gauge, level11, nGaugesTotal, nInflowGaugesTotal
    use mo_string_utils, only : num2str
    use mo_utils, only : ge

    implicit none

    character(256) :: fName

    integer(i4) :: i, iDomain, domainID, j

    integer(i4) :: err


    fName = trim(adjustl(dirConfigOut)) // trim(adjustl(file_config))
    call message()
    call message('  Log-file written to ', trim(fName))
    open(uconfig, file = fName, status = 'unknown', action = 'write', iostat = err)
    if (err .ne. 0) then
      call message('  Problems while creating File')
      call message('  Error-Code', num2str(err))
      stop
    end if
    write(uconfig, 200)
    write(uconfig, 100) 'mRM-UFZ v-' // trim(version)
    write(uconfig, 100) 'S. Thober, L. Samaniego & R. Kumar, UFZ'
    write(uconfig, 200)
    write(uconfig, 100)
    write(uconfig, 201) '         M A I N  mRM  C O N F I G U R A T I O N  I N F O R M A T I O N         '
    write(uconfig, 100)
    write(uconfig, 103) 'Number of domains            ', domainMeta%nDomains
    write(uconfig, 103) 'Total No. of gauges         ', nGaugesTotal
    write(uconfig, 103)    'Time Step [h]               ', timeStep
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      write(uconfig, 103) 'Total No. of nodes          ', level11(iDomain)%nCells
      write(uconfig, 103) 'No. of cells L0             ', level0(L0_Domain(iDomain))%nCells
      write(uconfig, 103) 'No. of cells L1             ', level1(iDomain)%nCells
      write(uconfig, 103) 'No. of cells L11            ', level11(iDomain)%nCells

      !    select case (iFlag_cordinate_sys)
      !    case (0)
      write(uconfig, 301)      'domain  ', domainID, '   Hydrology Resolution [m]      ', resolutionHydrology(iDomain)
      write(uconfig, 301)   'domain  ', domainID, '   Routing Resolution [m]        ', resolutionRouting(iDomain)
      !    case(1)
      !      write(uconfig, 302)       'domain  ',domainID, '   Hydrology Resolution [o]      ', resolutionHydrology(iDomain)
      !      write(uconfig, 302)   'domain  ',domainID, '   Routing Resolution [o]        ', resolutionRouting(iDomain)
      !    end select
    end do
    write(uconfig, 126)    'Flag READ  restart            ', read_restart
    write(uconfig, 126)    'Flag WRITE restart            ', write_restart
    !
    !******************
    ! Model Run period 
    !******************
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      write(uconfig, 115) '                      Model Run Periods for domain ', num2str(domainID)
      write(uconfig, 116) &
              'From                To', &
              '   Day Month  Year   Day Month  Year'
      write(uconfig, 117)  &
              'Warming Period (1)            ', &
              warmPer(iDomain)%dStart, warmPer(iDomain)%mStart, warmPer(iDomain)%yStart, &
              warmPer(iDomain)%dEnd, warmPer(iDomain)%mEnd, warmPer(iDomain)%yEnd, &
              'Evaluation Period (2)         ', &
              evalPer(iDomain)%dStart, evalPer(iDomain)%mStart, evalPer(iDomain)%yStart, &
              evalPer(iDomain)%dEnd, evalPer(iDomain)%mEnd, evalPer(iDomain)%yEnd, &
              'Simulation Period (1)+(2)     ', &
              SimPer(iDomain)%dStart, SimPer(iDomain)%mStart, SimPer(iDomain)%yStart, &
              SimPer(iDomain)%dEnd, SimPer(iDomain)%mEnd, SimPer(iDomain)%yEnd
    end do

    !*********************************
    ! Model Land Cover Observations 
    !*********************************
    if (processMatrix(8, 1) .eq. 1) then
      do iDomain = 1, domainMeta%nDomains
        domainID = domainMeta%indices(iDomain)
        write(uconfig, 118) '       Land Cover Observations for domain ', num2str(domainID)
        write(uconfig, 119) ' Start Year', ' End Year', '    Land cover scene', 'Land Cover File'
        do i = 1, nLCoverScene
          write(uconfig, 120) LC_year_start(i), LC_year_end(i), &
                  LCyearId(max(evalPer(iDomain)%yStart, LC_year_start(i)), iDomain), trim(LCfilename(i))
        end do
      end do
    end if
    !*********************************
    ! Initial Parameter Ranges
    !*********************************
    write(uconfig, 121) '  Initial Transfer Function Parameter Ranges (gammas)  '
    !
    ! Transfer functions
    write(uconfig, 122)      &
            '         i', '            min', '            max', '        current', &
            '                               name'
    do i = 1, size(global_parameters, 1)
      write(uconfig, 123) &
              i, global_parameters(i, 1), global_parameters(i, 2), global_parameters(i, 3), &
              trim(adjustl(global_parameters_name(i)))
    end do
    ! domain runoff data
    write(uconfig, 202) '                domain Runoff Data                '
    write(uconfig, 107) ' Gauge No.', '  domain Id', '     Qmax[m3/s]', '     Qmin[m3/s]'
    do i = 1, nGaugesTotal
      if(any(gauge%Q(:, i) > nodata_dp)) then
        write(uconfig, 108) i, gauge%domainId(i), maxval(gauge%Q(:, i), gauge%Q(:, i) > nodata_dp), &
                minval(gauge%Q(:, i), gauge%Q(:, i) > nodata_dp)
      else
        write(uconfig, 108) i, gauge%domainId(i), nodata_dp, nodata_dp
      end if
    end do
    ! inflow gauge data
    if (nInflowGaugesTotal .GT. 0) then
      write(uconfig, 202) '                domain Inflow Data                 '
      write(uconfig, 107) ' Gauge No.', '  domain Id', '     Qmax[m3/s]', '     Qmin[m3/s]'
      do i = 1, nInflowGaugesTotal
        if(all(InflowGauge%Q(:, i) > nodata_dp)) then
          write(uconfig, 108) i, InflowGauge%domainId(i), maxval(InflowGauge%Q(:, i), InflowGauge%Q(:, i) > nodata_dp), &
                  minval(InflowGauge%Q(:, i), InflowGauge%Q(:, i) > nodata_dp)
        else
          write(uconfig, 108) i, InflowGauge%domainId(i), nodata_dp, nodata_dp
        end if
      end do
    end if
    ! domain config
    write(uconfig, 218) 'domain-wise Configuration'
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      write(uconfig, 103) 'domain No.                   ', domainID, &
              'No. of gauges               ', basin_mrm(iDomain)%nGauges

      write(uconfig, 222)   'Directory list'

      write(uconfig, 224) 'Directory to morphological input         ', dirMorpho(iDomain)
      write(uconfig, 224) 'Directory to land cover input            ', dirLCover(iDomain)
      write(uconfig, 224) 'Directory to gauging station input       ', dirGauges(iDomain)
      if (mrm_coupling_mode .eq. 0) then
        write(uconfig, 224) 'Directory to simulated runoff input      ', dirTotalRunoff(iDomain)
      end if
      write(uconfig, 224) 'Directory to write output by default     ', dirOut(iDomain)
      write(uconfig, 224) 'Directory to write output when restarted ', dirRestartOut(iDomain)

      write(uconfig, 102) 'River Network  (Routing level)'
      write(uconfig, 100) 'Label 0 = intermediate draining cell '
      write(uconfig, 100) 'Label 1 = headwater cell             '
      write(uconfig, 100) 'Label 2 = sink cell                  '

      if (processMatrix(8, 1) .eq. 1_i4) then
        write(uconfig, 104) '   Overall', &
                '      From', &
                '        To', &
                '   Routing', &
                '     Label', &
                '    Length', &
                '      Mean', &
                '      Link', &
                '   Routing', &
                '   Routing', &
                '  Sequence', &
                '          ', &
                '          ', &
                '     Slope'
        !
        write(uconfig, 105) '        Id', &
                '      Node', &
                '      Node', &
                '', &
                '', &
                '      [km]', &
                '    [o/oo]'
        !
        do j = level11(iDomain)%iStart, level11(iDomain)%iEnd - 1
          i = L11_netPerm(j) + level11(iDomain)%iStart - 1 ! adjust permutation for multi-domain option
          write(uconfig, 106) i, L11_fromN(i), L11_toN(i), L11_rOrder(i), L11_label(i), &
                  L11_length(i) / 1000.0_dp, L11_slope(i) * 1.0e3_dp
        end do

      else if (processMatrix(8, 1) .eq. 2_i4) then
        write(uconfig, 134) '   Overall', &
                '      From', &
                '        To', &
                '   Routing', &
                '     Label', &
                '      Link', &
                '   Routing', &
                '   Routing', &
                '  Sequence', &
                '          '
        !
        write(uconfig, 135) '        Id', &
                '      Node', &
                '      Node', &
                '', &
                ''
        !
        do j = level11(iDomain)%iStart, level11(iDomain)%iEnd - 1
          i = L11_netPerm(j) + level11(iDomain)%iStart - 1 ! adjust permutation for multi-domain option
          write(uconfig, 136) i, L11_fromN(i), L11_toN(i), L11_rOrder(i), L11_label(i)
        end do
      end if
      ! draining node at L11
      write(uconfig, 109)  '   Overall', '     domain', &
              '      Cell', '   Routing', &
              '        Id', '   Node Id'
      do i = level11(iDomain)%Id(1), level11(iDomain)%Id(level11(iDomain)%nCells)
        write(uconfig, 110) i + level11(iDomain)%iStart - 1, i
      end do

      ! L1 level information
      write(uconfig, 111)  '  Modeling', '   Routing', ' Effective', &
              '      Cell', '   Cell Id', '      Area', &
              '        Id', '       [-]', '     [km2]'
      if (ge(resolutionRouting(iDomain), resolutionHydrology(iDomain))) then
        do i = level1(iDomain)%Id(1), level1(iDomain)%Id(level1(iDomain)%nCells)
          write(uconfig, 113) i + level1(iDomain)%iStart - 1, L1_L11_Id (i + level1(iDomain)%iStart - 1), &
            level1(iDomain)%CellArea(i) * 1.E-6_dp
      domainID = domainMeta%indices(iDomain)
        end do
      else
        do i = level11(iDomain)%Id(1), level11(iDomain)%Id(level11(iDomain)%nCells)
          write(uconfig, 110) i + level11(iDomain)%iStart - 1, L11_L1_Id (i + level11(iDomain)%iStart - 1)
        end do
      end if
      write(uconfig, 114)  ' Total[km2]', sum(level1(iDomain)%CellArea) * 1.E-6_dp
    end do

    write(uconfig, *)
    close(uconfig)

    !! Formats
    100 format (a80)
    102 format (/ 30('-') / a30 / 30('-'))
    103 format (a20, 10x, i10)
    104 format (/ 75('-') / 5a10, 5x, 2a10 / 5a10, 5x, 2a10)
    105 format (5a10, 5x, 2a10 / 75('-'))
    106 format (5i10, 5x, 2f10.3)
    107 format (2a10, 2a15)
    108 format (2i10, 2f15.3)
    !
    109 format (/ 20('-') / 2a10 / 2a10 / 2a10 / 20('-'))
    110 format (2i10)
    !
    111 format (/ 30('-') / 3a10 / 3a10 / 3a10 /  30('-'))
    113 format (            2i10,   1f10.3         )
    114 format (30('-') / a15, 5x, 1f10.3 /)
    !
    115 format (/61('-')/ a50, a10 /61('-'))
    116 format (39x, a22 / 25x, a36)
    117 format (3(a25, 6(i6)))
    !
    118 format (/50('-')/ a40, a10  /50('-'))
    119 format (a10, a10, a20, a20/)
    120 format (i10, i10, 10x, i10, a20)
    !
    121 format (/55('-')/ a55 /55('-'))
    122 format (a10, 3a15, a35)
    123 format (i10, 3f15.3, a35)
    !
    126 format (a30, 9x, L1)
    !
    134 format (/ 50('-') / 5a10 / 5a10)
    135 format (5a10 / 50('-'))
    136 format (5i10)
    !
    200 format (80('-'))
    201 format (a80)
    202 format (/50('-')/ a50 /50('-'))
    !
    218 format (/ 80('-')/ 26x, a24, 26x, /80('-'))
    222 format (/80('-')/ 26x, a21 /80('-'))
    224 format (a40, 5x, a256)

    301 format (a7, i2, a32, f15.0)
    ! 302 format (a7, i2, a32,es20.8)
  end Subroutine write_configfile

  ! ------------------------------------------------------------------

  !    NAME
  !        write_daily_obs_sim_discharge

  !    PURPOSE
  !>       \brief Write a file for the daily observed and simulated discharge timeseries
  !>       during the evaluation period for each gauging station

  !>       \details Write a file for the daily observed and simulated discharge timeseries
  !>       during the evaluation period for each gauging station

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:, :) :: Qobs" daily time series of observed dischargedims = (nModeling_days
  !>       , nGauges_total)
  !>       \param[in] "real(dp), dimension(:, :) :: Qsim" daily time series of modeled dischargedims = (nModeling_days ,
  !>       nGauges_total)

  !    HISTORY
  !>       \authors Rohini Kumar

  !>       \date August 2013

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine write_daily_obs_sim_discharge(Qobs, Qsim)

    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : dirOut, domainMeta
    use mo_errormeasures, only : kge, nse
    use mo_julian, only : dec2date
    use mo_message, only : message
    use mo_mrm_file, only : file_daily_discharge, ncfile_discharge, udaily_discharge
    use mo_mrm_global_variables, only : basin_mrm, gauge
    use mo_ncwrite, only : var2nc
    use mo_string_utils, only : num2str
    use mo_utils, only : ge

    implicit none

    ! daily time series of observed dischargedims = (nModeling_days , nGauges_total)
    real(dp), dimension(:, :), intent(in) :: Qobs

    ! daily time series of modeled dischargedims = (nModeling_days , nGauges_total)
    real(dp), dimension(:, :), intent(in) :: Qsim

    character(256) :: fName, formHeader, formData, dummy
    character(256), dimension(1) :: dnames

    integer(i4) :: domainID, iDomain, gg, tt, err

    integer(i4) :: igauge_start, igauge_end

    integer(i4) :: day, month, year

    integer(i4) :: tlength

    ! time axis
    integer(i4), allocatable, dimension(:) :: taxis

    real(dp) :: newTime

    logical :: create


    ! initalize igauge_start
    igauge_start = 1

    ! domain loop
    do iDomain = 1, domainMeta%nDomains
      domainID = domainMeta%indices(iDomain)
      if(basin_mrm(iDomain)%nGauges .lt. 1) cycle

      ! estimate igauge_end
      igauge_end = igauge_start + basin_mrm(iDomain)%nGauges - 1

      ! check the existance of file
      fName = trim(adjustl(dirOut(iDomain))) // trim(adjustl(file_daily_discharge))
      open(udaily_discharge, file = trim(fName), status = 'unknown', action = 'write', iostat = err)
      if(err .ne. 0) then
        call message ('  IOError while openening ', trim(fName))
        call message ('  Error-Code ', num2str(err))
        stop
      end if

      ! header
      write(formHeader, *) '( 4a8, ', basin_mrm(iDomain)%nGauges, '(2X, a5, i10.10, 2X, a5, i10.10) )'
      write(udaily_discharge, formHeader) 'No', 'Day', 'Mon', 'Year', &
              ('Qobs_', gauge%gaugeId(gg), &
                      'Qsim_', gauge%gaugeId(gg), gg = igauge_start, igauge_end)

      ! form data
      write(formData, *) '( 4I8, ', basin_mrm(iDomain)%nGauges, '(2X,   f15.7 , 2X,  f15.7  ) )'

      ! write data
      newTime = real(evalPer(iDomain)%julStart, dp) - 0.5_dp

      do tt = 1, (evalPer(iDomain)%julEnd - evalPer(iDomain)%julStart + 1)
        call dec2date(newTime, yy = year, mm = month, dd = day)
        write(udaily_discharge, formData) tt, day, month, year, (Qobs(tt, gg), Qsim(tt, gg), gg = igauge_start, igauge_end)
        newTime = newTime + 1.0_dp
      end do

      ! close file
      close(udaily_discharge)

      ! ======================================================================
      ! write netcdf file
      ! ======================================================================
      dnames(1) = 'time'
      ! dnames(2) = 'gauges'
      fName = trim(adjustl(dirOut(iDomain))) // trim(adjustl(ncfile_discharge))
      tlength = evalPer(iDomain)%julEnd - evalPer(iDomain)%julStart + 1
      create = .true.
      do gg = igauge_start, igauge_end
        ! write simulated discharge at that gauge
        call var2nc(trim(fName), Qsim(1 : tlength, gg), &
                dnames(1 : 1), 'Qsim_' // trim(num2str(gauge%gaugeID(gg), '(i10.10)')), create = create, &
                units = 'm3 s-1', long_name = 'simulated discharge at gauge ' // trim(num2str(gauge%gaugeID(gg), '(i10.10)')))
        create = .false.
        ! write observed discharge at that gauge
        call var2nc(trim(fName), Qobs(1 : tlength, gg), &
                dnames(1 : 1), 'Qobs_' // trim(num2str(gauge%gaugeID(gg), '(i10.10)')), create = create, &
                units = 'm3 s-1', long_name = 'observed discharge at gauge ' // trim(num2str(gauge%gaugeID(gg), '(i10.10)')))
      end do
      ! add time axis
      allocate(taxis(tlength))
      forall(tt = 1 : tlength) taxis(tt) = tt * 24 - 1
      call dec2date(real(evalPer(iDomain)%julStart, dp) - 0.5_dp, yy = year, mm = month, dd = day)
      call var2nc(trim(fName), taxis, &
              dnames(1 : 1), dnames(1), &
              units = 'hours since ' // &
                      trim(num2str(year)) // '-' // trim(num2str(month, '(i2.2)')) // '-' // trim(num2str(day, '(i2.2)')) // &
                      ' 00:00:00', &
              long_name = 'time in hours')
      deallocate(taxis)

      ! ======================================================================
      ! screen output
      ! ======================================================================
      call message()
      write(dummy, '(I3)') domainID
      call message('  OUTPUT: saved daily discharge file for domain ', trim(adjustl(dummy)))
      call message('    to ', trim(fname))
      do gg = igauge_start, igauge_end
        if (count(ge(Qobs(:, gg), 0.0_dp)) > 1)  then
          call message('    KGE of daily discharge (gauge #', trim(adjustl(num2str(gg))), '): ', &
                  trim(adjustl(num2str(kge(Qobs(:, gg), Qsim(:, gg), mask = (ge(Qobs(:, gg), 0.0_dp)))))))
          call message('    NSE of daily discharge (gauge #', trim(adjustl(num2str(gg))), '): ', &
                  trim(adjustl(num2str(nse(Qobs(:, gg), Qsim(:, gg), mask = (ge(Qobs(:, gg), 0.0_dp)))))))
        end if
      end do

      ! update igauge_start
      igauge_start = igauge_end + 1
      !
    end do
    !
  end subroutine write_daily_obs_sim_discharge

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_write_output_fluxes

  !    PURPOSE
  !>       \brief write fluxes to netcdf output files

  !>       \details This subroutine creates a netcdf data set
  !>       for writing L11_QTIN for different time averages.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"
  !>       \param[in] "integer(i4) :: nCells"
  !>       \param[in] "integer(i4) :: timeStep_model_outputs" timestep of model outputs
  !>       \param[in] "integer(i4) :: warmingDays"            number of warming days
  !>       \param[in] "real(dp) :: newTime"                   julian date of next time step
  !>       \param[in] "integer(i4) :: nTimeSteps"             number of total timesteps
  !>       \param[in] "integer(i4) :: nTStepDay"              number of timesteps per day
  !>       \param[in] "integer(i4) :: tt"                     current model timestep
  !>       \param[in] "integer(i4) :: day"                    current day of the year
  !>       \param[in] "integer(i4) :: month"                  current month of the year
  !>       \param[in] "integer(i4) :: year"                   current year
  !>       \param[in] "integer(i4) :: timestep"               current model time resolution
  !>       \param[in] "logical, dimension(:, :) :: mask11"    mask at level 11
  !>       \param[in] "real(dp), dimension(:) :: L11_qMod"    current routed streamflow

  !    HISTORY
  !>       \authors Stephan Thober

  !>       \date Aug 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_write_output_fluxes(iDomain, nCells, timeStep_model_outputs, warmingDays, newTime, nTimeSteps, &
                                    nTStepDay, tt, day, month, year, timestep, mask11, L11_qmod)

    use mo_julian, only : caldat
    use mo_kind, only : dp, i4

    implicit none

    integer(i4), intent(in) :: iDomain

    integer(i4), intent(in) :: nCells

    ! timestep of model outputs
    integer(i4), intent(in) :: timeStep_model_outputs

    ! number of warming days
    integer(i4), intent(in) :: warmingDays

    ! julian date of next time step
    real(dp), intent(in) :: newTime

    ! number of total timesteps
    integer(i4), intent(in) :: nTimeSteps

    ! number of timesteps per day
    integer(i4), intent(in) :: nTStepDay

    ! current model timestep
    integer(i4), intent(in) :: tt

    ! current day of the year
    integer(i4), intent(in) :: day

    ! current month of the year
    integer(i4), intent(in) :: month

    ! current year
    integer(i4), intent(in) :: year

    ! current model time resolution
    integer(i4), intent(in) :: timestep

    ! mask at level 11
    logical, intent(in), dimension(:, :), pointer :: mask11

    ! current routed streamflow
    real(dp), intent(in), dimension(:) :: L11_qMod

    integer(i4) :: tIndex_out

    logical :: writeout

    ! year of next timestep (newTime)
    integer(i4) :: new_year

    ! month of next timestep (newTime)
    integer(i4) :: new_month

    ! day of next timestep (newTime)
    integer(i4) :: new_day


    !
    if (tt .EQ. 1) then
      day_counter = day
      month_counter = month
      year_counter = year
      average_counter = 0_i4
    end if

    ! update the counters
    if (day_counter   .NE. day) day_counter = day
    if (month_counter .NE. month) month_counter = month
    if (year_counter  .NE. year)  year_counter = year
    call caldat(int(newTime), yy = new_year, mm = new_month, dd = new_day)

    ! output only for evaluation period
    tIndex_out = (tt - warmingDays * NTSTEPDAY) ! tt if write out of warming period

    if ((tIndex_out .gt. 0_i4)) then
      average_counter = average_counter + 1

      ! create output dataset
      if (tIndex_out .EQ. 1) nc = OutputDataset(iDomain, mask11, nCells)
      ! print*, 'After Init of OutputDatasetInit'

      ! update Dataset
      call nc%updateDataset(&
              1, &
              size(L11_Qmod), &
              L11_Qmod          &
              )

      ! determine write flag
      writeout = .false.
      if (timeStep_model_outputs .gt. 0) then
        if ((mod(tIndex_out, timeStep_model_outputs) .eq. 0) .or. (tt .eq. nTimeSteps)) writeout = .true.
      else
        select case(timeStep_model_outputs)
        case(0) ! only at last time step
          if (tt .eq. nTimeSteps) writeout = .true.
        case(-1) ! daily
          if ((day_counter .ne. new_day)     .or. (tt .eq. nTimeSteps)) writeout = .true.
        case(-2) ! monthly
          if ((month_counter .ne. new_month) .or. (tt .eq. nTimeSteps)) writeout = .true.
        case(-3) ! yearly
          if ((year_counter .ne. new_year)   .or. (tt .eq. nTimeSteps)) writeout = .true.
        case default ! no output at all

        end select
      end if

      ! write data
      if (writeout) call nc%writeTimestep(tIndex_out * timestep - 1)

      ! close dataset
      if (tt .eq. nTimeSteps) call nc%close()

    end if

  end subroutine mrm_write_output_fluxes

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_write_optifile

  !    PURPOSE
  !>       \brief Write briefly final optimization results.

  !>       \details Write overall best objective function and the best optimized parameter set to a file_opti.

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: best_OF"                             best objective function value as returnedby the
  !>       optimization routine
  !>       \param[in] "real(dp), dimension(:) :: best_paramSet"         best associated global parameter setCalled only
  !>       when optimize is .TRUE.
  !>       \param[in] "character(len = *), dimension(:) :: param_names"

  !    HISTORY
  !>       \authors David Schaefer

  !>       \date July 2013

  ! Modifications:
  ! Rohini Kumar   Aug 2013 - change in structure of the code including call statements
  ! Juliane Mai    Oct 2013 - clear parameter names added 
  !                         - double precision written
  ! Stephan Thober Oct 2015 - ported to mRM
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_write_optifile(best_OF, best_paramSet, param_names)

    use mo_common_mhm_mrm_file, only : file_opti, uopti
    use mo_common_variables, only : dirConfigOut
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    ! best objective function value as returnedby the optimization routine
    real(dp), intent(in) :: best_OF

    ! best associated global parameter setCalled only when optimize is .TRUE.
    real(dp), dimension(:), intent(in) :: best_paramSet

    character(len = *), dimension(:), intent(in) :: param_names

    character(256) :: fName, formHeader, formParams

    integer(i4) :: ii, err, n_params


    ! number of parameters
    n_params = size(best_paramSet)

    ! open file
    fName = trim(adjustl(dirConfigOut)) // trim(adjustl(file_opti))
    open(uopti, file = fName, status = 'unknown', action = 'write', iostat = err, recl = (n_params + 1) * 40)
    if(err .ne. 0) then
      call message ('  IOError while openening ', trim(fName))
      call message ('  Error-Code ', num2str(err))
      stop
    end if

    ! header 
    write(formHeader, *) '(a40,', n_params, 'a40)'
    ! len(param_names(1))=256 but only 39 characters taken here
    ! write(uopti, formHeader) 'OF', (trim(adjustl(param_names(ii))), ii=1, n_params)
    write(uopti, formHeader) 'OF', (trim(adjustl(param_names(ii)(1 : 39))), ii = 1, n_params)

    ! output
    write(formParams, *) '( es40.14, ', n_params, '(es40.14) )'
    write(uopti, formParams) best_OF, (best_paramSet(ii), ii = 1, n_params)

    ! close file
    close(uopti)

    ! screen output
    call message()
    call message(' Optimized parameters written to ', trim(fName))

  end subroutine mrm_write_optifile

  ! ------------------------------------------------------------------

  !    NAME
  !        mrm_write_optinamelist

  !    PURPOSE
  !>       \brief Write final, optimized parameter set in a namelist format.

  !>       \details Write final, optimized parameter set in a namelist format.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:, :) :: parameters"                               (min, max, opti)
  !>       \param[in] "logical, dimension(size(parameters, 1)) :: maskpara"                   .true. if parameter was
  !>       calibrated
  !>       \param[in] "character(len = *), dimension(size(parameters, 1)) :: parameters_name" clear names of parameters

  !    HISTORY
  !>       \authors Juliane Mai

  !>       \date Dec 2013

  ! Modifications:
  ! Stephan Thober Oct 2015 - adapted to mRM
  ! Stephan Thober Nov 2016 - adapt header to routing process
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine mrm_write_optinamelist(parameters, maskpara, parameters_name)

    use mo_common_mhm_mrm_file, only : file_opti_nml, uopti_nml
    use mo_common_variables, only : dirConfigOut, processMatrix
    use mo_message, only : message
    use mo_string_utils, only : num2str

    implicit none

    ! (min, max, opti)
    real(dp), dimension(:, :), intent(in) :: parameters

    ! .true. if parameter was calibrated
    logical, dimension(size(parameters, 1)), intent(in) :: maskpara

    ! clear names of parameters
    character(len = *), dimension(size(parameters, 1)), intent(in) :: parameters_name

    character(256) :: fName

    character(3) :: flag

    integer(i4) :: err

    integer(i4) :: iPar


    ! open file
    fName = trim(adjustl(dirConfigOut)) // trim(adjustl(file_opti_nml))
    open(uopti_nml, file = fName, status = 'unknown', action = 'write', iostat = err)
    if(err .ne. 0) then
      call message ('  IOError while openening ', trim(fName))
      call message ('  Error-Code ', num2str(err))
      stop
    end if

    write(uopti_nml, *) '!global_parameters'
    write(uopti_nml, *) '!PARAMETER                       lower_bound  upper_bound          value   FLAG  SCALING'

    write(uopti_nml, *) '! ', trim(adjustl('routing'))

    if (processMatrix(8, 1) .eq. 1_i4) write(uopti_nml, *) '&routing1'
    if (ProcessMatrix(8, 1) .eq. 2_i4) write(uopti_nml, *) '&routing2'
    if (ProcessMatrix(8, 1) .eq. 3_i4) write(uopti_nml, *) '&routing3'

    do iPar = 1, size(parameters, 1)
      if (maskpara(iPar)) then
        flag = ' 1 '
      else
        flag = ' 0 '
      end if
      write(uopti_nml, *) trim(adjustl(parameters_name(iPar))), ' = ', &
              parameters(iPar, 1), ' , ', &
              parameters(iPar, 2), ' , ', &
              parameters(iPar, 3), ' , ', &
              flag, ', 1 '
    end do

    write(uopti_nml, *) '/'
    write(uopti_nml, *) ' '

    ! close file
    close(uopti_nml)

    ! screen output
    call message()
    call message(' Optimized parameters written in namelist format to ', trim(fName))

  end subroutine mrm_write_optinamelist


end module mo_mrm_write
