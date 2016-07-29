!> \file mo_meteo_forcings.f90

!> \brief Prepare meteorological forcings data for mHM.

!> \details Prepare meteorological forcings data for mHM.

!> \authors Rohini Kumar
!> \date Jan 2012

MODULE mo_meteo_forcings

  ! This module provides routines to read meteorological data.

  ! Written  Rohini Kumar, Jan 2013
  ! Modified

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prepare_meteo_forcings_data  
 
  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         prepare_meteo_forcings_data
  
  !     PURPOSE
  !>        \brief Prepare meteorological forcings data for a given variable

  !>        \details Prepare meteorological forcings data for a given variable.
  !>                 Internally this subroutine calls another routine meteo_wrapper   
  !>                 for different meterological variables

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in] "integer(i4)              :: iBasin"        Basin Id

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

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date Jan 2013
  !           Matthias Zink,   Jun  2013 - addded NetCDf reader
  !           Rohini Kumar,    Aug  2013 - name changed "inputFormat" to inputFormat_meteo_forcings
  !           Matthias Zink,   Feb  2014 - added read in for different PET processes (process 5)
  !           Stephan Thober,  Jun  2014 - add chunk_config for chunk read, 
  !                                        copied L2 initialization to mo_startup
  !
  subroutine prepare_meteo_forcings_data(iBasin, tt)
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_timer,            only:                         &
         timer_start, timer_stop, timer_get                   ! Timing of processes
    use mo_global_variables, only: &
         dirPrecipitation, dirTemperature,                  & ! directory of meteo input
         dirReferenceET,                                    & ! PET input path  if process 5 is 'PET is input' (case 0)
         dirMinTemperature, dirMaxTemperature,              & ! PET input paths if process 5 is Hargreaves-Samani (case 1)
         dirNetRadiation,                                   & ! PET input paths if process 5 is Priestley-Taylor (case 2)
         dirabsVapPressure, dirwindspeed,                   & ! PET input paths if process 5 is Penman-Monteith (case 3)
         inputFormat_meteo_forcings,                        & ! 'bin' for binary data or 'nc' for NetCDF input
         nBasins,                                           & ! Number of basins for multi-basin optimization 
         processMatrix,                                     & ! process configuration
         readPer, timeStep_model_inputs,                    & ! chunk read in config                           
         L1_pre, L1_temp, L1_pet , L1_tmin, L1_tmax,        & ! meteorological data
         L1_netrad, L1_absvappress, L1_windspeed              ! meteorological data

    implicit none

    integer(i4),                   intent(in)  :: iBasin        ! Basin Id
    integer(i4),                   intent(in)  :: tt            ! current timestep

    ! local variables
    logical                                    :: read_flag     ! indicate whether data should be read

    ! configuration of chunk_read
    call chunk_config( iBasin, tt, read_flag, readPer )
 
    ! only read, if read_flag is true
    if ( read_flag ) then
       ! free L1 variables if chunk read is activated
       if ( timeStep_model_inputs(iBasin) .ne. 0 ) then
          if ( allocated( L1_pre         )) deallocate( L1_pre         )
          if ( allocated( L1_temp        )) deallocate( L1_temp        )
          if ( allocated( L1_pet         )) deallocate( L1_pet         )
          if ( allocated( L1_tmin        )) deallocate( L1_tmin        )
          if ( allocated( L1_tmax        )) deallocate( L1_tmax        )
          if ( allocated( L1_netrad      )) deallocate( L1_netrad      )
          if ( allocated( L1_absvappress )) deallocate( L1_absvappress )
          if ( allocated( L1_windspeed   )) deallocate( L1_windspeed   )
       end if
       
       !  basin characteristics and read meteo header
       if ( timeStep_model_inputs(iBasin) .eq. 0 ) then
          call message( '  Reading meteorological forcings for basin: ', trim(adjustl(num2str(iBasin))),' ...')
          call timer_start(1)
       end if

       ! precipitation
       if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read precipitation        ...' )
       call meteo_forcings_wrapper( iBasin, dirPrecipitation(iBasin), inputFormat_meteo_forcings, &
            L1_pre, lower=0.0_dp, upper=1000._dp, ncvarName='pre' )

       ! temperature
       if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read temperature          ...' )
       call meteo_forcings_wrapper( iBasin, dirTemperature(iBasin), inputFormat_meteo_forcings,  &
            L1_temp, lower = -100._dp, upper=100._dp, ncvarName='tavg' )

       ! read input for PET (process 5) depending on specified option
       ! 0 - input, 1 - Hargreaves-Samani, 2 - Priestley-Taylor, 3 - Penman-Monteith
       select case (processMatrix(5,1))    

       case(0) ! pet is input
          if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read pet                  ...' )
          call meteo_forcings_wrapper( iBasin, dirReferenceET(iBasin), inputFormat_meteo_forcings, &
               L1_pet, lower=0.0_dp, upper = 1000._dp, ncvarName='pet' )
          ! allocate PET and dummies for mhm_call
          if ((iBasin.eq.nBasins) .OR. (timeStep_model_inputs(iBasin) .NE. 0)) then
             allocate( L1_tmin(1,1)); allocate( L1_tmax(1,1) ); allocate( L1_netrad(1,1) )
             allocate( L1_absvappress(1,1)); allocate( L1_windspeed(1,1) )
          end if

       case(1) ! Hargreaves-Samani formulation (input: minimum and maximum Temperature)
          if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read min. temperature     ...' )
          call meteo_forcings_wrapper( iBasin, dirMinTemperature(iBasin), inputFormat_meteo_forcings, &
               L1_tmin, lower=-100.0_dp, upper = 100._dp, ncvarName='tmin' )
          if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read max. temperature     ...' )
          call meteo_forcings_wrapper( iBasin, dirMaxTemperature(iBasin), inputFormat_meteo_forcings, &
               L1_tmax, lower=-100.0_dp, upper = 100._dp, ncvarName='tmax' )
          ! allocate PET and dummies for mhm_call
          if ((iBasin .eq. nBasins) .OR. (timeStep_model_inputs(iBasin) .NE. 0)) then
             allocate( L1_pet    (size(L1_tmax, dim=1), size(L1_tmax, dim=2)))
             allocate( L1_netrad(1,1) ); allocate( L1_absvappress(1,1)); allocate( L1_windspeed(1,1) )
          end if

       case(2) ! Priestley-Taylor formulation (input: net radiation)
          if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read net radiation        ...' )
          call meteo_forcings_wrapper( iBasin, dirNetRadiation(iBasin), inputFormat_meteo_forcings, &
               L1_netrad, lower=-500.0_dp, upper = 1500._dp, ncvarName='net_rad' )
          ! allocate PET and dummies for mhm_call
          if ((iBasin .eq. nBasins) .OR. (timeStep_model_inputs(iBasin) .NE. 0)) then
             allocate( L1_pet    (size(L1_netrad, dim=1), size(L1_netrad, dim=2)))
             allocate( L1_tmin(1,1)); allocate( L1_tmax(1,1) )
             allocate( L1_absvappress(1,1)); allocate( L1_windspeed(1,1) )
          end if

       case(3) ! Penman-Monteith formulation (input: net radiationm absulute vapour pressure, windspeed)
          if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read net radiation        ...' )
          call meteo_forcings_wrapper( iBasin, dirNetRadiation(iBasin), inputFormat_meteo_forcings, &
               L1_netrad, lower=-500.0_dp, upper = 1500._dp, ncvarName='net_rad' )
          if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read absolute vapour pressure  ...' )
          call meteo_forcings_wrapper( iBasin, dirabsVapPressure(iBasin), inputFormat_meteo_forcings, &
               L1_absvappress, lower=0.0_dp, upper = 15000.0_dp, ncvarName='eabs' )
          if ( timeStep_model_inputs(iBasin) .eq. 0 ) call message( '    read windspeed            ...' )
          call meteo_forcings_wrapper( iBasin, dirwindspeed(iBasin), inputFormat_meteo_forcings, &
               L1_windspeed, lower=0.0_dp, upper = 250.0_dp, ncvarName='windspeed' )
          ! allocate PET and dummies for mhm_call
          if ((iBasin.eq.nBasins) .OR. (timeStep_model_inputs(iBasin) .NE. 0)) then
             allocate( L1_pet    (size(L1_absvappress, dim=1), size(L1_absvappress, dim=2)))
             allocate( L1_tmin(1,1)); allocate( L1_tmax(1,1) )
          end if
       end select

       if ( timeStep_model_inputs(iBasin) .eq. 0 ) then
          call timer_stop(1)
          call message('    in ', trim(num2str(timer_get(1),'(F9.3)')), ' seconds.')
       end if
    end if

end subroutine prepare_meteo_forcings_data



  ! ------------------------------------------------------------------

  !     NAME
  !         meteo_forcings_wrapper
  
  !     PURPOSE
  !>        \brief Prepare meteorological forcings data for mHM at Level-1

  !>        \details Prepare meteorological forcings data for mHM, which include \n
  !>         1) Reading meteo. datasets at their native resolution for every basin \n
  !>         2) Perform aggregation or disaggregation of meteo. datasets from their \n
  !>            native resolution (level-2) to the required hydrologic resolution (level-1)\n
  !>         3) Pad the above datasets of every basin to their respective global ones
  !>                 

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in] "integer(i4)               :: iBasin"        Basin Id
  !>        \param[in] "character(len=*)          :: dataPath"      Data path where a given meteo. variable is stored

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[in] "real(dp), dimension(:,:)  :: dataOut1"      Packed meterological variable for the whole simulation period

  !     INTENT(IN), OPTIONAL
  !>        \param[in] "real(dp), optional        :: lower"    Lower bound for check of validity of data values
  !>        \param[in] "real(dp), optional        :: upper"    Upper bound for check of validity of data values
  !>        \param[in] "type(period), optional    :: readPer"  reading Period
  !>        \param[in] "character(len=*), optional:: ncvarName" name of the variable (for .nc files)


  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !         None

  !     RESTRICTIONS

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Rohini Kumar
  !>        \date Jan 2013
  !         Modified, Stephan Thober, Jun 2014 -- changed to readPer
  !                   Stephan Thober, Feb 2016 -- refactored deallocate statements

  subroutine meteo_forcings_wrapper(iBasin, dataPath, inputFormat, dataOut1, lower, upper, ncvarName)
  
    use mo_global_variables,           only: readPer, level1, level2
    use mo_mhm_constants,              only: nodata_dp
    use mo_init_states,                only: get_basin_info
    use mo_read_meteo,                 only: read_meteo_bin
    use mo_read_forcing_nc,            only: read_forcing_nc
    use mo_spatial_agg_disagg_forcing, only: spatial_aggregation, spatial_disaggregation
    use mo_append,                     only: append                    ! append vector
    
    implicit none

    integer(i4),                          intent(in)    :: iBasin        ! Basin Id
    character(len=*),                     intent(in)    :: dataPath      ! Data path
    character(len=*),                     intent(in)    :: inputFormat   ! either 'bin' or 'nc'
    real(dp), dimension(:,:),allocatable, intent(inout) :: dataOut1      ! Packed meteorological variable
    real(dp),                   optional, intent(in)    :: lower         ! lower bound for data points
    real(dp),                   optional, intent(in)    :: upper         ! upper bound for data points
    character(len=*),           optional, intent(in)    :: ncvarName     ! name of the variable (for .nc files)

    
    integer(i4)                                :: nrows1, ncols1
    logical, dimension(:,:), allocatable       :: mask1
    integer(i4)                                :: ncells1

    integer(i4)                                :: nrows2, ncols2
    logical, dimension(:,:), allocatable       :: mask2

    real(dp), dimension(:,:,:), allocatable    :: L2_data            ! meteo data at level-2 
    real(dp), dimension(:,:,:), allocatable    :: L1_data            ! meteo data at level-1
    real(dp), dimension(:,:), allocatable      :: L1_data_packed     ! packed meteo data at level-1 from 3D to 2D

    integer(i4)                                :: nTimeSteps
    real(dp)                                   :: cellFactorHbyM   ! level-1_resolution/level-2_resolution
    integer(i4)                                :: t

    ! get basic basin information at level-1
    call get_basin_info( iBasin, 1, nrows1, ncols1, nCells=nCells1, mask=mask1 )
    
    ! make  basic basin information at level-2
    call get_basin_info( iBasin, 2, nrows2, ncols2, mask=mask2 )

    select case (trim(inputFormat))       
    case ('bin')
       ! read data
       if( present(lower) .AND. (.not. present(upper)) ) then
          CALL read_meteo_bin( dataPath, nRows2, nCols2, readPer,  L2_data, mask2, lower=lower )
       end if
       !
       if( present(upper) .AND. (.not. present(lower)) ) then
          CALL read_meteo_bin( dataPath, nRows2, nCols2, readPer, L2_data, mask2, upper=upper )
       end if
       !
       if( present(lower) .AND. present(upper) ) then
          CALL read_meteo_bin( dataPath, nRows2, nCols2, readPer, L2_data, mask2, lower=lower, upper=upper )
       end if
    
       if( (.not. present(lower)) .AND. (.not. present(upper)) ) then
          CALL read_meteo_bin( dataPath, nRows2, nCols2, readPer, L2_data, mask2 )
       end if
    case('nc')
       if( present(lower) .AND. (.not. present(upper)) ) then
          CALL read_forcing_nc( dataPath, nRows2, nCols2, readPer, ncvarName, L2_data, mask2, &
               lower=lower )
       end if
       !
       if( present(upper) .AND. (.not. present(lower)) ) then
          CALL read_forcing_nc( dataPath, nRows2, nCols2, readPer, ncvarName, L2_data, mask2, &
               upper=upper )
       end if
       !
       if( present(lower) .AND. present(upper) ) then
          CALL read_forcing_nc( dataPath, nRows2, nCols2, readPer, ncvarName, L2_data, mask2, &
               lower=lower, upper=upper )
       end if
    
       if( (.not. present(lower)) .AND. (.not. present(upper)) ) then
          CALL read_forcing_nc( dataPath, nRows2, nCols2, readPer, ncvarName, L2_data, mask2 )
       end if
    case DEFAULT
       stop '***ERROR: meteo_forcings_wrapper: Not recognized input format'
    end select

    ! cellfactor to decide on the upscaling or downscaling of meteo. fields
    cellFactorHbyM = level1%cellsize(iBasin) / level2%cellsize(iBasin) 

    ! upscaling & packing
    if(cellFactorHbyM .gt. 1.0_dp) then 
        call spatial_aggregation(L2_data, level2%cellsize(iBasin), level1%cellsize(iBasin), mask1, mask2, L1_data)
    ! downscaling   
    elseif(cellFactorHbyM .lt. 1.0_dp) then
        call spatial_disaggregation(L2_data, level2%cellsize(iBasin), level1%cellsize(iBasin), mask1, mask2, L1_data)
    ! nothing
    else
      allocate( L1_data( size(L2_data,1), size(L2_data,2), size(L2_data,3) ) )
      L1_data(:,:,:) = L2_data(:,:,:)
   end if
   ! free memory immediately
   deallocate(L2_data)
    
    ! pack variables
    nTimeSteps = size(L1_data, 3)
    allocate( L1_data_packed(nCells1, nTimeSteps))
    do t = 1, nTimeSteps
       L1_data_packed(:,t) = pack( L1_data(:,:,t), MASK=mask1(:,:) ) 
    end do
    ! free memory immediately
    deallocate(L1_data)
    
    ! append
    call append( dataOut1, L1_data_packed(:,:), nodata_dp )

    !free space
    deallocate(L1_data_packed) 
    
  end subroutine meteo_forcings_wrapper

  ! ------------------------------------------------------------------
  !
  ! subroutine chunk_config
  !
  ! determines the start date, end date, and read_flag 
  ! given basin id and current timestep
  !
  ! author: Stephan Thober
  !
  ! created: June 2014
  ! ------------------------------------------------------------------
  subroutine chunk_config( iBasin, tt, read_flag, readPer )
    !
    use mo_kind,             only: i4
    use mo_global_variables, only: period
    use mo_mhm_constants,    only: nodata_dp
    !
    implicit none
    !
    ! input variables
    integer(i4), intent(in)  :: iBasin ! current Basin
    integer(i4), intent(in)  :: tt     ! current timestep
    !
    ! output variables
    logical,     intent(out) :: read_flag  ! indicate whether reading data should be read
    type(period),intent(out) :: readPer    ! start and end dates of reading Period

    ! initialize
    read_flag        = .false.
    if ( tt .eq. 1_i4 ) then
       readPer%julStart = int( nodata_dp, i4 )
       readPer%julEnd   = int( nodata_dp, i4 )
       readPer%dstart   = int( nodata_dp, i4 )
       readPer%mstart   = int( nodata_dp, i4 )
       readPer%ystart   = int( nodata_dp, i4 )
       readPer%dend     = int( nodata_dp, i4 )
       readper%mend     = int( nodata_dp, i4 )
       readper%yend     = int( nodata_dp, i4 )
       readPer%Nobs     = int( nodata_dp, i4 )
    end if

    ! evaluate date and timeStep_model_inputs to get read_flag -------
    read_flag = is_read( iBasin, tt )
    !
    ! determine start and end date of chunk to read
    if ( read_flag ) call chunk_size( iBasin, tt, readPer )
    !
  end subroutine chunk_config
  ! ------------------------------------------------------------------

  !     NAME
  !         is_read
  
  !     PURPOSE
  !>        \brief evaluate whether new chunk should be read at this timestep

  !     CALLING SEQUENCE
  !         flag = is_read( iBasin, tt )

  !     INTENT(IN)
  !>        \param[in] "integer(i4)              :: iBasin"    current Basin
  !>        \param[in] "integer(i4)              :: tt"        current time step

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

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Jun 2014
  ! ------------------------------------------------------------------
  function is_read( iBasin, tt )
    
    use mo_kind,             only: i4
    use mo_global_variables, only: simPer, timeStep_model_inputs, timestep, nTstepDay
    use mo_message,          only: message
    use mo_julian,           only: caldat

    ! input variables
    integer(i4), intent(in) :: iBasin ! Basin ID
    integer(i4), intent(in) :: tt     ! timestep
    
    ! return variable
    logical :: is_read

    ! local variables
    integer(i4)             :: Ndays        ! number of simulated days
    integer(i4)             :: day          ! day
    integer(i4)             :: month        ! months
    integer(i4)             :: year         ! years
    integer(i4)             :: Ndays_before ! number of simulated days one timestep before
    integer(i4)             :: day_before   ! day one simulated timestep before
    integer(i4)             :: month_before ! month one simulated timestep before
    integer(i4)             :: year_before  ! year one simulated timestep before
    
    
    ! initialize
    is_read      = .false.
    
    ! special case for first timestep
    if ( tt .eq. 1_i4 ) then
       is_read = .true.
    else
       ! check if a new day started by comparing the day of the current time step (Ndays)
       ! with the one before (Ndays_before)
       Ndays        = ceiling( (real(tt,dp)           ) / real(nTstepDay, dp) ) 
       Ndays_before = ceiling( (real(tt,dp) - 1.0_dp  ) / real(nTstepDay, dp) )   
       
       ! evaluate cases of given timeStep_model_inputs
       select case( timeStep_model_inputs(iBasin) )
       case(0)  ! only at the beginning of the period
          if ( tt .eq. 1_i4 ) is_read = .true.
       case(1:) ! every timestep with frequency timeStep_model_inputs
          if ( mod( (tt-1) * timestep, timeStep_model_inputs(iBasin) * 24 ) .eq. 0_i4 ) is_read = .true.
       case(-1) ! every day
          if ( Ndays .ne. Ndays_before ) is_read = .true.
       case(-2) ! every month
          if ( Ndays .ne. Ndays_before ) then
             ! calculate months
             call caldat( simPer(iBasin)%julStart + Ndays - 1, dd = day, mm = month, yy = year )
             call caldat( simPer(iBasin)%julStart + Ndays_before - 1, dd = day_before, mm = month_before, yy = year_before )
             if ( month .ne. month_before ) is_read = .true.
          end if
       case(-3) ! every year
          if ( Ndays .ne. Ndays_before ) then
             ! calculate months
             call caldat( simPer(iBasin)%julStart + Ndays - 1, dd = day, mm = month, yy = year )
             call caldat( simPer(iBasin)%julStart + Ndays_before - 1, dd = day_before, mm = month_before, yy = year_before )
             if ( year .ne. year_before ) is_read = .true.
          end if
       case default ! not specified correctly
          call message('ERROR*** mo_meteo_forcings: function is_read: timStep_model_inputs not specified correctly!')
          stop
       end select
    end if
  
  end function is_read
  ! ------------------------------------------------------------------

  !     NAME
  !         chunk_size
  
  !     PURPOSE
  !>        \brief calculate beginning and end of read Period, i.e. that
  !>               is length of current chunk to read

  !     CALLING SEQUENCE
  !         call chunk_size( iBasin, tt, readPer )

  !     INTENT(IN)
  !>        \param[in] "integer(i4)              :: iBasin"    current Basin to process
  !>        \param[in] "integer(i4)              :: tt"        current time step
  !>        \param[in] "type(period)             :: readPer"   start and end dates of read Period

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

  !     EXAMPLE

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Jun 2014
  !         modified Stephan Thober - Jan 2015 added iBasin
  ! ------------------------------------------------------------------
  subroutine chunk_size( iBasin, tt, readPer )
    
    use mo_kind,             only: i4
    use mo_global_variables, only: simPer, & !     start and end of simulation period
         timeStep_model_inputs,            & !     frequency for reading meteo input  
         period,                           & !     type for describing read in period
         nTstepDay                           !     Number of time intervals per day
    use mo_message,          only: message
    use mo_julian,           only: caldat, julday
    
    implicit none
    
    ! input variables
    integer(i4), intent(in)  :: tt         ! current time step
    integer(i4), intent(in)  :: iBasin     ! Basin ID
    
    ! output variables
    type(period),intent(out) :: readPer    ! start and end dates of reading Period

    ! local variables
    integer(i4)              :: Ndays        ! number of simulated days
    integer(i4)              :: day          ! day
    integer(i4)              :: month        ! months
    integer(i4)              :: year         ! years

    ! calculate date of start date
    Ndays        =  ceiling( real(tt,dp) / real(nTstepDay,dp) )

    ! get start date
    readPer%julStart = simPer(iBasin)%julStart + Ndays - 1
    
    ! calculate end date according to specified frequency
    select case ( timeStep_model_inputs(iBasin) )
    case(0)  ! length of chunk has to cover whole period
       readPer%julEnd = simPer(iBasin)%julEnd 
    case(1:) ! every timestep with frequency timeStep_model_inputs
       readPer%julEnd= readPer%julStart + timeStep_model_inputs(iBasin) - 1
    case(-1) ! every day
       readPer%julEnd = readPer%julStart
    case(-2) ! every month
       ! calculate date
       call caldat( simPer(iBasin)%julStart + Ndays, dd = day, mm = month, yy = year )
       ! increment month
       if ( month .eq. 12 ) then
          month = 1
          year  = year + 1
       else
          month = month + 1
       end if
       readPer%julEnd = julday( dd = 1, mm = month, yy = year ) - 1
    case(-3) ! every year
       ! calculate date
       call caldat( simPer(iBasin)%julStart + Ndays, dd = day, mm = month, yy = year )
       readPer%julEnd = julday( dd = 31, mm = 12, yy = year )
    case default ! not specified correctly
       call message('ERROR*** mo_meteo_forcings: chunk_size: timStep_model_inputs not specified correctly!')
       stop
    end select

    ! end date should not be greater than end of simulation period
    readPer%julEnd = min( readPer%julEnd, simPer(iBasin)%julEnd )
    
    ! calculate the dates of the start and end dates
    call caldat( readPer%julStart, dd = readPer%dstart, mm = readPer%mstart, yy = readPer%ystart )
    call caldat( readPer%julEnd,   dd = readPer%dEnd,   mm = readPer%mend,   yy = readPer%yend   )
    readPer%Nobs = readPer%julEnd - readPer%julstart + 1    
    
  end subroutine chunk_size
  !
END MODULE mo_meteo_forcings
