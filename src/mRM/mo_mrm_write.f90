!> \file mo_mrm_write.f90

!> \brief write of discharge and restart files

!> \details This module contains the subroutines for
!> writing the discharge files and optionally the restart
!> files.

!> \author Stephan Thober
!> \date Aug 2015
module mo_mrm_write

  use mo_kind, only: i4, dp
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

  !     NAME
  !         mrm_write

  !     PURPOSE
  !>        \brief write discharge and restart files
  !
  !>        \details First, this subroutine calls the writing or restart files that only
  !>        succeeds if it happens after the write of mHM restart files because
  !>        mHM restart files must exist. Second, simulated discharge is aggregated to the daily
  !>        scale and then written to file jointly with observed discharge
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
  !         None
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
  !         None
  !
  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Juliane Mai, Rohini Kumar & Stephan Thober
  !>        \date Aug 2015
  !         Modified, 

  subroutine mrm_write()

    use mo_mrm_global_variables, only: &
         mRM_runoff, &
         gauge, nGaugesTotal, basin_mrm, nBasins, evalPer, warmingDays_mrm, simPer, &
         ntstepday, write_restart, dirRestartOut, mrm_coupling_mode
    use mo_mrm_restart,  only: mrm_write_restart

    implicit none
    
    ! local variables
    integer(i4) :: iBasin
    integer(i4) :: iDay, iS, iE
    integer(i4) :: ii
    integer(i4) :: tt
    integer(i4) :: gg
    integer(i4) :: nTimeSteps
    real(dp), dimension(:,:), allocatable :: d_Qmod

    ! --------------------------------------------------------------------------
    ! WRITE CONFIG FILE
    ! --------------------------------------------------------------------------
    if (mrm_coupling_mode .ne. 2) call write_configfile()

    ! --------------------------------------------------------------------------
    ! WRITE RESTART
    ! --------------------------------------------------------------------------
    if (write_restart) then
       do iBasin = 1, nBasins
          call mrm_write_restart(iBasin, dirRestartOut)
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
    ii = maxval( evalPer(1:nBasins)%julEnd - evalPer(1:nBasins)%julStart + 1 )
    allocate( d_Qmod(ii, nGaugesTotal) ) 
    d_Qmod = 0.0_dp

    ! loop over basins
    do ii = 1, nBasins
       nTimeSteps = ( simPer(ii)%julEnd - simPer(ii)%julStart + 1 ) * NTSTEPDAY
       iDay = 0
       ! loop over timesteps
       do tt = warmingDays_mrm(ii)*NTSTEPDAY+1, nTimeSteps, NTSTEPDAY
          iS = tt
          iE = tt + NTSTEPDAY - 1
          iDay = iDay + 1
          ! over gauges
          do gg = 1, basin_mrm%nGauges(ii)
             d_Qmod(iDay, basin_mrm%gaugeIndexList(ii,gg) ) = &
                  sum( mRM_runoff(iS:iE, basin_mrm%gaugeIndexList(ii,gg)) )/ real(NTSTEPDAY,dp)
          end do
          !
       end do
    end do
    ! write in an ASCII file          ! OBS[nModeling_days X nGauges_total] , SIM[nModeling_days X nGauges_total] 
    call write_daily_obs_sim_discharge( gauge%Q(:,:), d_Qmod(:,:) )
    ! free space
    deallocate(d_Qmod)        
  end subroutine mrm_write
  ! ------------------------------------------------------------------

  !      NAME
  !         write_configfile

  !     PURPOSE
  !>        \brief This modules writes the results of the configuration into an ASCII-file

  !>        \details 

  !     INTENT(IN)
  !         None

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
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         

  !     LITERATURE
  !         

  !     HISTORY
  !>        \author Christoph Schneider
  !>        \date May 2013
  !         Modified, Juliane Mai,    May 2013 - module version and documentation
  !                   Stephan Thober, Jun 2014 - bug fix in L11 config print out 
  !                   Stephan Thober, Jun 2014 - updated read_restart
  !                   Rohini, Luis  , Jul 2015 - updated version, L1 level prints
  !                   Stephan Thober, Sep 2015 - updated write of stream network

    Subroutine write_configfile()

    use mo_utils,            only: ge
    use mo_kind,             only: i4, dp
    use mo_mrm_constants,    only: nodata_dp
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    USE mo_mrm_file,         only: file_config, uconfig, version
    use mo_mrm_global_variables, only: &
         mrm_coupling_mode,          &
         nBasins,                    &
         basin_mrm,                  &
         gauge,                      &
         InflowGauge,                &
         L11_nCells,                 &
         L11_netPerm,                &
         L11_fromN,                  &
         L11_toN,                    &
         L11_rOrder,                 &
         L11_label,                  &
         L11_length,                 &         
         L11_slope,                  &
         L11_ID,                     &
         L11_L1_Id,                  &
         L1_areaCell,                &
         L1_L11_Id,                  &
         L0_nCells,                  &
         L1_nCells,                  &
         nGaugesTotal,               &
         nInflowGaugesTotal,         &
         resolutionRouting,          &  
         dirGauges,                  &
         dirTotalRunoff,             &
         timeStep,                   &
         resolutionHydrology,        &
         read_restart,               &
         write_restart,              &
         dirConfigOut,               &
         dirMorpho,                  &
         dirLCover,                  &
         dirOut,                     &
         dirRestartOut,              &  
         warmPer,                    &
         evalPer,                    &
         SimPer,                     &
         LCyearId,                   &
         LCfilename

    use mo_common_variables, only: &
         global_parameters,      &
         global_parameters_name 

    implicit none
    !
    ! local
    !
    character(256)                         :: fName
    integer(i4)                            :: i, j, n
    integer(i4)                            :: err

    fName=  trim(adjustl(dirConfigOut))//trim(adjustl(file_config))
    call message()
    call message('  Log-file written to ', trim(fName))
    open(uconfig, file=fName, status='unknown', action='write', iostat=err)
    if (err .ne. 0) then
       call message('  Problems while creating File' )
       call message('  Error-Code', num2str(err) )
       stop
    end if
    write(uconfig, 200) 
    write(uconfig, 100) 'mRM-UFZ v-'//trim(version)
    write(uconfig, 100) 'S. Thober, L. Samaniego & R. Kumar, UFZ'
    write(uconfig, 200) 
    write(uconfig, 100)
    write(uconfig, 201) '         M A I N  mRM  C O N F I G U R A T I O N  I N F O R M A T I O N         '
    write(uconfig, 100)
    write(uconfig, 103) 'Number of basins            ', nBasins
    write(uconfig, 103) 'Total No. of nodes          ', L11_nCells
    write(uconfig, 103) 'Total No. of reaches        ', L11_nCells-1
    write(uconfig, 103) 'No. of cells L0             ', L0_nCells
    write(uconfig, 103) 'No. of cells L1             ', L1_nCells
    write(uconfig, 103) 'No. of cells L11            ', L11_nCells
    write(uconfig, 103) 'Total No. of gauges         ', nGaugesTotal
    write(uconfig, 103)    'Time Step [h]               ', timeStep
    do i=1, nBasins
    !    select case (iFlag_cordinate_sys)
    !    case (0)
          write(uconfig, 301)      'Basin  ',i, '   Hydrology Resolution [m]      ', resolutionHydrology(i)
          write(uconfig, 301)   'Basin  ',i, '   Routing Resolution [m]        ', resolutionRouting(i)
    !    case(1)
    !      write(uconfig, 302)       'Basin  ',i, '   Hydrology Resolution [o]      ', resolutionHydrology(i)
    !      write(uconfig, 302)   'Basin  ',i, '   Routing Resolution [o]        ', resolutionRouting(i)
    !    end select
    end do
    write(uconfig, 126)    'Flag READ  restart            ', read_restart
    write(uconfig, 126)    'Flag WRITE restart            ', write_restart
    !
    !******************
    ! Model Run period 
    !******************
    do j = 1, nBasins
       write(uconfig, 115) '                      Model Run Periods for Basin ', num2str(j)
       write(uconfig, 116) &
            'From                To', &
            '   Day Month  Year   Day Month  Year'
       write(uconfig,117)  &
            'Warming Period (1)            ',&
            warmPer(j)%dStart, warmPer(j)%mStart, warmPer(j)%yStart  ,& 
            warmPer(j)%dEnd  , warmPer(j)%mEnd  , warmPer(j)%yEnd    ,&
            'Evaluation Period (2)         ',&
            evalPer(j)%dStart   ,evalPer(j)%mStart   , evalPer(j)%yStart      ,& 
            evalPer(j)%dEnd     ,evalPer(j)%mEnd     , evalPer(j)%yEnd        ,&
            'Simulation Period (1)+(2)     ',&
            SimPer(j)%dStart  , SimPer(j)%mStart  , SimPer(j)%yStart  ,&
            SimPer(j)%dEnd    , SimPer(j)%mEnd    , SimPer(j)%yEnd
    end do

    !*********************************
    ! Model Land Cover Observations 
    !*********************************
    do j = 1, nBasins
       write(uconfig,118) '       Land Cover Observations for Basin ', num2str(i)
       write(uconfig,119) '      Year', '    Land cover scene', 'Land Cover File'
       do i=1,SimPer(j)%yEnd-SimPer(j)%yStart+1
          write(uconfig,120) i+SimPer(j)%yStart-1, LCyearId(i+SimPer(j)%yStart-1, j), &
               trim(LCfilename(LCyearId(i+SimPer(j)%yStart-1, j)))
       end do
    end do
    !*********************************
    ! Initial Parameter Ranges
    !*********************************
    write(uconfig,121) '  Initial Transfer Function Parameter Ranges (gammas)  '
    !
    ! Transfer functions
    write(uconfig,122)      &
         '         i', '            min', '            max', '        current', &
         '                               name'
    do i=1, size(global_parameters,1)
       write(uconfig,123) &
            i, global_parameters(i,1), global_parameters(i,2), global_parameters(i,3), &
            trim(adjustl(global_parameters_name(i)))
    end do
    ! basin runoff data
    write(uconfig, 202) '                Basin Runoff Data                '
    write(uconfig, 107) ' Gauge No.', '  Basin Id', '     Qmax[m3/s]', '     Qmin[m3/s]'         
    do i=1, nGaugesTotal
       if( any(gauge%Q(:,i) > nodata_dp) ) then
          write(uconfig,108) i, gauge%basinId(i), maxval(gauge%Q(:,i), gauge%Q(:,i) > nodata_dp), &
               minval(gauge%Q(:,i), gauge%Q(:,i) > nodata_dp)
       else
          write(uconfig,108) i, gauge%basinId(i), nodata_dp, nodata_dp
       end if
    end do
    ! inflow gauge data
    if ( nInflowGaugesTotal .GT. 0 ) then
       write(uconfig, 202) '                Basin Inflow Data                 '
       write(uconfig, 107) ' Gauge No.', '  Basin Id', '     Qmax[m3/s]', '     Qmin[m3/s]'         
       do i=1, nInflowGaugesTotal
          if( all(InflowGauge%Q(:,i) > nodata_dp) ) then
             write(uconfig,108) i, InflowGauge%basinId(i), maxval(InflowGauge%Q(:,i), InflowGauge%Q(:,i) > nodata_dp), &
                  minval(InflowGauge%Q(:,i), InflowGauge%Q(:,i) > nodata_dp)
          else
             write(uconfig,108) i, InflowGauge%basinId(i), nodata_dp, nodata_dp
          end if
       end do
    end if
    ! basin config
    write(uconfig,218) 'Basin-wise Configuration'
    do n=1,nBasins
       write(uconfig,103) 'Basin No.                   ', n, &
            'No. of gauges               ', basin_mrm%nGauges(n)

       write(uconfig, 222)   'Directory list'

       write(uconfig, 224) 'Directory to morphological input         ',  dirMorpho(n)
       write(uconfig, 224) 'Directory to land cover input            ',  dirLCover(n)
       write(uconfig, 224) 'Directory to gauging station input       ',  dirGauges(n)
       if (mrm_coupling_mode .eq. 0) then
          write(uconfig, 224) 'Directory to simulated runoff input      ',  dirTotalRunoff(n)
       end if
       write(uconfig, 224) 'Directory to write output by default     ',  dirOut(n)
       write(uconfig, 224) 'Directory to write output when restarted ',  dirRestartOut(n)

       write(uconfig, 102) 'River Network  (Routing level)'
       write(uconfig, 100) 'Label 0 = intermediate draining cell '
       write(uconfig, 100) 'Label 1 = headwater cell             '
       write(uconfig, 100) 'Label 2 = sink cell                  '

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
            '',           &
            '      [km]', &
            '    [o/oo]'
       !
       do j=basin_mrm%L11_iStart(n), basin_mrm%L11_iEnd(n)-1
          i=L11_netPerm(j) + basin_mrm%L11_iStart(n) - 1 ! adjust permutation for multi-basin option
          write(uconfig,106) i, L11_fromN(i), L11_toN(i), L11_rOrder(i), L11_label(i), &
               L11_length(i)/1000.0_dp, L11_slope(i)*1.0e3_dp
       end do
       ! draining node at L11
       write(uconfig, 109)  '   Overall', '     Basin', &
            '      Cell', '   Routing', &
            '        Id', '   Node Id'
       do i=basin_mrm%L11_iStart(n), basin_mrm%L11_iEnd(n)
          write(uconfig, 110) i, L11_Id(i)
       end do
       
       ! L1 level information
       write(uconfig, 111)  '  Modeling', '   Routing', ' Effective', &
            '      Cell', '   Cell Id', '      Area', &
            '        Id', '       [-]', '     [km2]'
       if (ge(resolutionRouting(n), resolutionHydrology(n))) then
          do i=basin_mrm%L1_iStart(n), basin_mrm%L1_iEnd(n)
             write(uconfig,113) i, L1_L11_Id(i), L1_areaCell(i)
          end do
       else
          do i=basin_mrm%L11_iStart(n), basin_mrm%L11_iEnd(n)
             write(uconfig,113) i, L11_L1_Id(i) !, L1_areaCell(i)
          end do
       end if
       write(uconfig,114)  ' Total[km2]', sum(L1_areaCell(basin_mrm%L1_iStart(n): basin_mrm%L1_iEnd(n)))
    end do

    write(uconfig,*)
    close(uconfig)  

    !! Formats
100 format (a80)
102 format (/ 30('-') / a30 / 30('-') )
103 format (a20, 10x, i10)
104 format (/ 75('-') / 5a10, 5x, 2a10 / 5a10, 5x, 2a10)
105 format (5a10, 5x, 2a10 / 75('-'))
106 format (5i10, 5x, 2f10.3 )
107 format (2a10, 2a15)
108 format (2i10, 2f15.3)
    !
109 format (/ 20('-') / 2a10 / 2a10 / 2a10 / 20('-') )
110 format (            2i10 )
    !
111 format (/ 30('-') / 3a10 / 3a10 / 3a10 /  30('-') )
113 format (            2i10,   1f10.3         )
114 format (30('-') / a15, 5x,  1f10.3 /       )
    !
115 format (/61('-')/ a50, a10 /61('-'))
116 format (39x,a22 / 25x, a36)
117 format ( 3(a25,6(i6) /) )
    !
118 format (/50('-')/ a40, a10  /50('-'))
119 format (a10,      a20, a20/)
120 format (i10, 10x, i10, a20)
    !
121 format (/55('-')/ a55 /55('-'))
122 format (a10, 3a15,   a35)
123 format (i10, 3f15.3, a35)
    !
126 format (a30,9x,L1)
    !
200 format (80('-'))
201 format (a80)
202 format (/50('-')/ a50 /50('-'))
    !
218 format (/ 80('-')/ 26x, a24,26x,  /80('-'))
222 format (/80('-')/ 26x,a21 /80('-'))
224 format (a40, 5x, a256)

301 format (a7, i2, a32,f15.0)
! 302 format (a7, i2, a32,es20.8)
  end Subroutine write_configfile

  ! ------------------------------------------------------------------

  !     NAME
  !         write_daily_obs_sim_discharge

  !     PURPOSE
  !>        \brief Write a file for the daily observed and simulated discharge timeseries 
  !>                during the evaluation period for each gauging station

  !>        \details Write a file for the daily observed and simulated discharge timeseries 
  !>                during the evaluation period for each gauging station

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in]  "real(dp), dimension(:,:)    :: Qobs"    daily time series of observed discharge 
  !>                                                             dims = (nModeling_days , nGauges_total)
  !>        \param[in]  "real(dp), dimension(:,:)    :: Qsim"    daily time series of modeled discharge
  !>                                                             dims = (nModeling_days , nGauges_total)

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
  !         None 

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date August 2013

  subroutine write_daily_obs_sim_discharge(Qobs, Qsim)


    use mo_errormeasures, only: kge, nse
    use mo_julian, only: dec2date
    use mo_message, only: message
    use mo_string_utils, only: num2str
    use mo_utils, only: ge
    use mo_mrm_file, only: file_daily_discharge, udaily_discharge, &
         ncfile_discharge
    use mo_mrm_global_variables, only: &
         nBasins, &
         basin_mrm, &
         dirOut, evalPer, &
         gauge
    use mo_ncwrite, only: var2nc

    implicit none

    ! input arguments
    real(dp), dimension(:,:), intent(in) :: Qobs      ! observed time series  [nModeling_days X nGauges_total]
    real(dp), dimension(:,:), intent(in) :: Qsim      ! simulated time series [nModeling_days X nGauges_total]

    ! local vars
    character(256) :: fName, formHeader, formData, dummy, dnames(1)
    integer(i4) :: bb, gg, tt, err
    integer(i4) :: igauge_start, igauge_end
    integer(i4) :: day, month, year
    integer(i4) :: tlength
    integer(i4), allocatable :: taxis(:) ! time axis
    real(dp) :: newTime
    logical :: create


    ! initalize igauge_start
    igauge_start = 1

    ! basin loop
    do bb = 1, nBasins
       if( basin_mrm%nGauges(bb) .lt. 1 ) cycle

       ! estimate igauge_end
       igauge_end = igauge_start + basin_mrm%nGauges(bb) - 1

       ! check the existance of file
       fName = trim(adjustl(dirOut(bb))) // trim(adjustl(file_daily_discharge))
       open(udaily_discharge, file=trim(fName), status='unknown', action='write', iostat=err)
       if( err .ne. 0 ) then
          call message ('  IOError while openening ',trim(fName))
          call message ('  Error-Code ', num2str(err))
          stop
       end if

       ! header
       write(formHeader, *) '( 4a8, ' , basin_mrm%nGauges(bb),'(2X, a5, i10.10, 2X, a5, i10.10) )' 
       write(udaily_discharge, formHeader) 'No', 'Day', 'Mon', 'Year', &
            ( 'Qobs_', gauge%gaugeId(gg), &
            'Qsim_', gauge%gaugeId(gg), gg=igauge_start, igauge_end )

       ! form data
       write(formData, *) '( 4I8, ' , basin_mrm%nGauges(bb),'(2X,   f15.7 , 2X,  f15.7  ) )' 

       ! write data
       newTime  = real(evalPer(bb)%julStart,dp) - 0.5_dp

       do tt = 1, (evalPer(bb)%julEnd - evalPer(bb)%julStart + 1)          
          call dec2date(newTime, yy=year, mm=month, dd=day)
          write(udaily_discharge, formData) tt, day, month, year, ( Qobs(tt,gg), Qsim(tt,gg) , gg=igauge_start, igauge_end )
          newTime = newTime + 1.0_dp
       end do

       ! close file
       close(udaily_discharge)

       ! ======================================================================
       ! write netcdf file
       ! ======================================================================
       dnames(1) = 'time'
       ! dnames(2) = 'gauges'
       fName = trim(adjustl(dirOut(bb))) // trim(adjustl(ncfile_discharge))
       tlength = evalPer(bb)%julEnd - evalPer(bb)%julStart + 1
       create = .true.
       do gg = igauge_start, igauge_end
          ! write simulated discharge at that gauge
          call var2nc(trim(fName), Qsim(1:tlength, gg), &
               dnames(1:1), 'Qsim_' // trim(num2str(gauge%gaugeID(gg), '(i5.5)')), create=create, &
               units='m3 s-1', long_name='simulated discharge at gauge ' // trim(num2str(gauge%gaugeID(gg), '(i5.5)')))
          create = .false.
          ! write observed discharge at that gauge
          call var2nc(trim(fName), Qobs(1:tlength, gg), &
               dnames(1:1), 'Qobs_' // trim(num2str(gauge%gaugeID(gg), '(i5.5)')), create=create, &
               units='m3 s-1', long_name='observed discharge at gauge ' // trim(num2str(gauge%gaugeID(gg), '(i5.5)')))
       end do
       ! add time axis
       allocate(taxis(tlength))
       forall(tt = 1:tlength) taxis(tt) = tt * 24 - 1
       call dec2date(real(evalPer(bb)%julStart,dp) - 0.5_dp, yy=year, mm=month, dd=day)
       call var2nc(trim(fName), taxis, &
            dnames(1:1), dnames(1), &
            units='hours since '// &
            trim(num2str(year))//'-'//trim(num2str(month, '(i2.2)'))//'-'//trim(num2str(day, '(i2.2)'))// &
            ' 00:00:00', &
            long_name='time in hours')
       deallocate(taxis)

       ! ======================================================================
       ! screen output
       ! ======================================================================
       call message()
       write(dummy,'(I3)') bb
       call message('  OUTPUT: saved daily discharge file for basin ', trim(adjustl(dummy)))
       call message('    to ',trim(fname))
       do gg=igauge_start, igauge_end
          call message('    KGE of daily discharge (gauge #',trim(adjustl(num2str(gg))),'): ', &
               trim(adjustl(num2str(kge(Qobs(:,gg), Qsim(:,gg), mask=(ge(Qobs(:,gg), 0.0_dp)))))) )
          call message('    NSE of daily discharge (gauge #',trim(adjustl(num2str(gg))),'): ', &
               trim(adjustl(num2str(nse(Qobs(:,gg), Qsim(:,gg), mask=(ge(Qobs(:,gg), 0.0_dp)))))) )
       end do

       ! update igauge_start
       igauge_start = igauge_end + 1
       !
    end do
    !
  end subroutine write_daily_obs_sim_discharge

  ! ------------------------------------------------------------------

  !     NAME
  !         mrm_write_output_fluxes

  !     PURPOSE
  !>        \brief write fluxes to netcdf output files
  !
  !>        \details This subroutine creates a netcdf data set
  !>        for writing L11_QTIN for different time averages.
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
  !         None
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
  !         None
  !
  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Aug 2015
  !         Modified, 

  subroutine mrm_write_output_fluxes( &
       iBasin, &
       timeStep_model_outputs, & ! timestep of model outputs
       warmingDays_mrm, & ! number of warming days
       newTime, & ! julian date of next time step
       nTimeSteps, & ! number of total timesteps
       nTStepDay, & ! number of timesteps per day
       tt, & ! current model timestep
       day, & ! current day of the year
       month, & ! current month of the year
       year, & ! current year
       timestep, & ! current model time resolution
       mask11, & ! mask at level 11
       L11_qmod & ! current routed streamflow
       )
    use mo_kind, only: i4, dp
    use mo_julian, only: caldat
    implicit none
    ! input variables
    integer(i4), intent(in) :: iBasin
    integer(i4), intent(in) :: timeStep_model_outputs
    integer(i4), intent(in) :: warmingDays_mrm
    real(dp), intent(in) :: newTime
    integer(i4), intent(in) :: nTimeSteps
    integer(i4), intent(in) :: nTStepDay
    integer(i4), intent(in) :: tt
    integer(i4), intent(in) :: day
    integer(i4), intent(in) :: month
    integer(i4), intent(in) :: year
    integer(i4), intent(in) :: timestep
    logical, intent(in) :: mask11(:,:)
    real(dp), intent(in) :: L11_qMod(:)
    ! local variables
    integer(i4) :: tIndex_out
    logical :: writeout
    integer(i4) :: new_year ! year of next timestep (newTime)
    integer(i4) :: new_month ! month of next timestep (newTime)
    integer(i4) :: new_day ! day of next timestep (newTime)
    !
    if ( tt .EQ. 1 ) then
       day_counter   = day
       month_counter = month
       year_counter  = year
       average_counter = 0_i4
    end if
    
    ! update the counters
    if (day_counter   .NE. day  ) day_counter   = day
    if (month_counter .NE. month) month_counter = month
    if (year_counter  .NE. year)  year_counter  = year
    call caldat(int(newTime), yy=new_year, mm=new_month, dd=new_day)

    ! output only for evaluation period
    tIndex_out = (tt-warmingDays_mrm*NTSTEPDAY) ! tt if write out of warming period

    if ((tIndex_out .gt. 0_i4)) then
       average_counter = average_counter + 1

       ! create output dataset 
       if ( tIndex_out .EQ. 1 ) nc = OutputDataset(iBasin, mask11)
       
       ! update Dataset
       call nc%updateDataset( &
            1               , &
            size(L11_Qmod)  , &
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
             if (((tIndex_out .gt. 1) .and. (day_counter .ne. new_day)) .or. (tt .eq. nTimeSteps))     writeout = .true.
          case(-2) ! monthly
             if (((tIndex_out .gt. 1) .and. (month_counter .ne. new_month)) .or. (tt .eq. nTimeSteps)) writeout = .true.
          case(-3) ! yearly
             if (((tIndex_out .gt. 1) .and. (year_counter .ne. new_year)) .or. (tt .eq. nTimeSteps))   writeout = .true.
          case default ! no output at all
             continue
          end select
       endif

       ! write data
       if (writeout) call nc%writeTimestep(tIndex_out*timestep-1)

       ! close dataset
       if (tt .eq. nTimeSteps) call nc%close()

    end if
    
  end subroutine mrm_write_output_fluxes

  ! ------------------------------------------------------------------

  !     NAME
  !         write_optifile

  !     PURPOSE
  !>        \brief Write briefly final optimization results.

  !>        \details Write overall best objective function and the best optimized parameter set to a file_opti.

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in]  "real(dp)                   :: best_OF"         best objective function value as returned 
  !>                                                                    by the optimization routine
  !>        \param[in]  "real(dp), dimension(:)     :: best_paramSet"   best associated global parameter set

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
  !         None

  !     RESTRICTIONS
  !>        Called only when optimize is .TRUE.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author David Schaefer
  !>        \date July 2013

  !         Modified, Rohini Kumar,   Aug 2013 - change in structure of the code including call statements
  !                   Juliane Mai,    Oct 2013 - clear parameter names added
  !                                            - double precision written
  !                   Stephan Thober, Oct 2015 - ported to mRM

  subroutine mrm_write_optifile(best_OF, best_paramSet, param_names)

    use mo_mrm_file,             only: file_opti, uopti
    use mo_mrm_global_variables, only: dirConfigOut
    use mo_message,              only: message
    use mo_string_utils,         only: num2str

    implicit none

    real(dp),                       intent(in)     :: best_OF     
    real(dp),         dimension(:), intent(in)     :: best_paramSet
    character(len=*), dimension(:), intent(in)     :: param_names

    ! local variables
    character(256)                         :: fName, formHeader, formParams
    integer(i4)                            :: ii, err, n_params

    ! number of parameters
    n_params = size(best_paramSet)

    ! open file
    fName = trim(adjustl(dirConfigOut)) // trim(adjustl(file_opti))
    open(uopti, file=fName, status='unknown', action='write', iostat=err, recl=(n_params+1)*40)
    if( err .ne. 0 ) then
       call message ('  IOError while openening ',trim(fName))
       call message ('  Error-Code ', num2str(err))
       stop
    end if

    ! header 
    write(formHeader, *) '(a40,',n_params,'a40)'
    ! len(param_names(1))=256 but only 39 characters taken here
    ! write(uopti, formHeader) 'OF', (trim(adjustl(param_names(ii))), ii=1, n_params)
    write(uopti, formHeader) 'OF', (trim(adjustl(param_names(ii)(1:39))), ii=1, n_params)

    ! output
    write(formParams, *) '( es40.14, ', n_params,'(es40.14) )' 
    write(uopti, formParams) best_OF, (best_paramSet(ii), ii=1, n_params)

    ! close file
    close(uopti)

    ! screen output
    call message()
    call message(' Optimized parameters written to ', trim(fName) )

  end subroutine mrm_write_optifile

  ! ------------------------------------------------------------------

  !     NAME
  !         write_optinamelist

  !     PURPOSE
  !>        \brief Write final, optimized parameter set in a namelist format.

  !>        \details  Write final, optimized parameter set in a namelist format. 

  !     CALLING SEQUENCE
  !         None

  !     INTENT(IN)
  !>        \param[in]  "real(dp)         :: parameters(:,:)"      information about parameter (min, max, opti)
  !>        \param[in]  "logical          :: maskpara(:)"          infomation which parameter where optimized
  !>        \param[in]  "character(len=*) :: parameters_name(:)"   clear names of parameters

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
  !         None

  !     RESTRICTIONS
  !>        Called only when optimize is .TRUE.

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2013

  !         Modified,
  !               Oct 2015, Stephan Thober - adapted to mRM

  subroutine mrm_write_optinamelist(parameters, maskpara, parameters_name)

    use mo_mrm_file,             only: file_opti_nml, uopti_nml
    use mo_mrm_global_variables, only: dirConfigOut
    use mo_message,              only: message
    use mo_string_utils,         only: num2str

    implicit none

    real(dp),         dimension(:,:),                intent(in) :: parameters        ! (min, max, opti)
    logical,          dimension(size(parameters,1)), intent(in) :: maskpara          ! .true. if parameter was calibrated
    character(len=*), dimension(size(parameters,1)), intent(in) :: parameters_name   ! clear names of parameters

    ! local variables
    character(256)                           :: fName
    character(3)                             :: flag
    integer(i4)                              :: err
    integer(i4)                              :: iPar

    ! open file
    fName = trim(adjustl(dirConfigOut)) // trim(adjustl(file_opti_nml))
    open(uopti_nml, file=fName, status='unknown', action='write', iostat=err)
    if( err .ne. 0 ) then
       call message ('  IOError while openening ',trim(fName))
       call message ('  Error-Code ', num2str(err))
       stop
    end if

    write(uopti_nml,*) '!global_parameters'
    write(uopti_nml,*) '!PARAMETER                       lower_bound  upper_bound          value   FLAG  SCALING'

    write(uopti_nml,*) '! ', trim(adjustl('routing'))

    write(uopti_nml,*) '&routing1'
    
    do iPar=1, size(parameters,1)
       if (maskpara(iPar)) then
          flag=' 1 '
       else
          flag=' 0 '
       end if
       write(uopti_nml,*) trim(adjustl(parameters_name(iPar))), ' = ', &
            parameters(iPar,1), ' , ', &
            parameters(iPar,2), ' , ', &
            parameters(iPar,3), ' , ', &
            flag, ', 1 '
    end do

    write(uopti_nml,*) '/'
    write(uopti_nml,*) ' '

    ! close file
    close(uopti_nml)

    ! screen output
    call message()
    call message(' Optimized parameters written in namelist format to ', trim(fName) )

  end subroutine mrm_write_optinamelist
    

end module mo_mrm_write
