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
  !           Rohini Kumar,   Aug  2013 - name changed "inputFormat" to inputFormat_meteo_forcings
  !
  subroutine prepare_meteo_forcings_data(iBasin, tt)
    use mo_message,          only: message
    use mo_string_utils,     only: num2str
    use mo_timer,            only:                         &
         timer_start, timer_stop, timer_get              ! Timing of processes
    use mo_global_variables, only: &
         dirPrecipitation, dirTemperature, dirReferenceET, &
         inputFormat_meteo_forcings,                       &
         L1_pre, L1_temp, L1_pet 

    implicit none

    integer(i4),                   intent(in)  :: iBasin        ! Basin Id
    integer(i4),                   intent(in)  :: tt            ! current timestep
    
    ! local variables
    logical                                    :: read_flag     ! indicate whether data should be read
    integer(i4)                                :: start_date    ! julian date of beginning of read period
    integer(i4)                                :: end_date      ! julian date of end of read period

    ! configuration of chunk_read
    call chunk_config(iBasin, tt, read_flag, start_date, end_date )
    read_flag = .FALSe.
    if ( tt .eq. 1 ) read_flag = .true.

    ! only read, if read_flag is true
    if ( read_flag ) then
       ! basic basin characteristics and read meteo header
       call message( '  Reading meteorological forcings for basin: ', trim(adjustl(num2str(iBasin))),' ...')
       call timer_start(1)

       call L2_variable_init(iBasin)

       ! precipitation
       call message( '    read precipitation ...' )
       call meteo_forcings_wrapper( iBasin, dirPrecipitation(iBasin), inputFormat_meteo_forcings, &
            L1_pre, lower=0.0_dp, upper=1000._dp, ncvarName='pre' )

       ! temperature
       call message( '    read temperature   ...' )
       call meteo_forcings_wrapper( iBasin, dirTemperature(iBasin), inputFormat_meteo_forcings,  &
            L1_temp, lower = -100._dp, upper=100._dp, ncvarName='tavg' )
    
       ! pet
       call message( '    read pet           ...' )
       call meteo_forcings_wrapper( iBasin, dirReferenceET(iBasin), inputFormat_meteo_forcings, &
            L1_pet, lower=0.0_dp, upper = 1000._dp, ncvarName='pet' )
       
       call timer_stop(1)
       call message('    in ', trim(num2str(timer_get(1),'(F9.3)')), ' seconds.')
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
  !>        \param[in] "real(dp), optional)           :: lower"    Lower bound for check of validity of data values
  !>        \param[in] "real(dp), optional)           :: upper"    Upper bound for check of validity of data values


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

  subroutine meteo_forcings_wrapper(iBasin, dataPath, inputFormat, dataOut1, lower, upper, ncvarName)
  
    use mo_global_variables,           only: simPer, level1, level2
    use mo_init_states,                only: get_basin_info
    use mo_read_meteo,                 only: read_meteo_bin, read_meteo_nc
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
          CALL read_meteo_bin( dataPath, nRows2, nCols2, simPer,  L2_data, mask2, lower=lower)
       end if
       !
       if( present(upper) .AND. (.not. present(lower)) ) then
          CALL read_meteo_bin( dataPath, nRows2, nCols2, simPer, L2_data, mask2, upper=upper)
       end if
       !
       if( present(lower) .AND. present(upper) ) then
          CALL read_meteo_bin( dataPath, nRows2, nCols2, simPer, L2_data, mask2, lower=lower, upper=upper)
       end if
    
       if( (.not. present(lower)) .AND. (.not. present(upper)) ) then
          CALL read_meteo_bin( dataPath, nRows2, nCols2, simPer, L2_data, mask2)
       end if
    case('nc')
       if( present(lower) .AND. (.not. present(upper)) ) then
          CALL read_meteo_nc( dataPath, nRows2, nCols2, simPer, ncvarName, L2_data, mask2, lower=lower)
       end if
       !
       if( present(upper) .AND. (.not. present(lower)) ) then
          CALL read_meteo_nc( dataPath, nRows2, nCols2, simPer, ncvarName, L2_data, mask2, upper=upper)
       end if
       !
       if( present(lower) .AND. present(upper) ) then
          CALL read_meteo_nc( dataPath, nRows2, nCols2, simPer, ncvarName, L2_data, mask2, lower=lower, upper=upper)
       end if
    
       if( (.not. present(lower)) .AND. (.not. present(upper)) ) then
          CALL read_meteo_nc( dataPath, nRows2, nCols2, simPer, ncvarName, L2_data, mask2)
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
    
    ! pack variables
    nTimeSteps = size(L1_data, 3)
    allocate( L1_data_packed(nCells1, nTimeSteps))
    do t = 1, nTimeSteps
       L1_data_packed(:,t) = pack( L1_data(:,:,t), MASK=mask1(:,:) ) 
    end do
    
    ! append
    call append( dataOut1, L1_data_packed(:,:) )

    !free space
    deallocate(L1_data, L2_data, L1_data_packed) 
    
  end subroutine meteo_forcings_wrapper

  ! ------------------------------------------------------------------

  !     NAME
  !         L2_variable_init
  
  !     PURPOSE
  !>        \brief Initalize Level-2 meteorological forcings data

  !>        \details following tasks are performed
  !>                 1)  cell id & numbering
  !>                 2)  mask creation
  !>                 3)  append variable of intrest to global ones

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
  !>        \date Feb 2013

  ! --------------------------------------------------------------------------
  subroutine L2_variable_init(iBasin)

    use mo_read_spatial_data,only: read_header_ascii
    use mo_message,          only: message
    use mo_append,           only: append                      ! append vector
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_init_states,      only: calculate_grid_properties

    use mo_global_variables, only: nBasins, basin, level0, level2, dirPrecipitation
    use mo_mhm_constants,    only: nodata_dp
    use mo_file,             only: file_meteo_header, umeteo_header

    
    implicit none

    integer(i4), intent(in)                   :: iBasin

    ! local variables
    integer(i4)                               :: nrows0, ncols0
    logical, dimension(:,:), allocatable      :: mask0
    real(dp)                                  :: xllcorner0, yllcorner0
    real(dp)                                  :: cellsize0

    integer(i4)                               :: nrows2, ncols2
    logical, dimension(:,:), allocatable      :: mask2
    real(dp)                                  :: xllcorner2, yllcorner2
    real(dp)                                  :: cellsize2
    integer(i4)                               :: nCells2
    real(dp)                                  :: cellFactor
    integer(i4)                               :: i, j, ic, jc
    character(256)                            :: fName

    !--------------------------------------------------------
    ! STEPS::
    ! 1) Estimate each variable locally for a given basin
    ! 2) Pad each variable to its corresponding global one
    !--------------------------------------------------------
    
    ! assign space
    if(iBasin .eq. 1) then
       allocate( level2%nrows        (nBasins) )
       allocate( level2%ncols        (nBasins) )
       allocate( level2%xllcorner    (nBasins) )
       allocate( level2%yllcorner    (nBasins) )
       allocate( level2%cellsize     (nBasins) )
       allocate( level2%nodata_value (nBasins) )
     end if

    ! read header file 
    ! NOTE: assuming the header file for all metero variables are same as that of precip.
    !       A counter check for this assumption is perfromed in the read_meteo_bin file 
    
    fName =  trim(adjustl(dirPrecipitation(iBasin))) // trim(adjustl(file_meteo_header))
    call read_header_ascii( trim(fName), umeteo_header,   &
                            level2%nrows(iBasin), level2%ncols(iBasin), level2%xllcorner(iBasin), &
                            level2%yllcorner(iBasin), level2%cellsize(iBasin),  level2%nodata_value(iBasin) )
  
   ! level-0 information
   call get_basin_info( iBasin, 0, nrows0, ncols0, mask=mask0,                         &
                        xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0 ) 
   ! grid information
   call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp,          &
                                   level2%cellsize(iBasin), &
                                   nrows2, ncols2, xllcorner2, yllcorner2, cellsize2,level2%nodata_value(iBasin) )

   ! check
   if (  (ncols2     .ne.  level2%ncols(iBasin))         .or. &
         (nrows2     .ne.  level2%nrows(iBasin))         .or. &
         ( abs(xllcorner2 - level2%xllcorner(iBasin)) .gt. tiny(1.0_dp) )     .or. &
         ( abs(yllcorner2 - level2%yllcorner(iBasin)) .gt. tiny(1.0_dp) )     .or. &
         ( abs(cellsize2  - level2%cellsize(iBasin))  .gt. tiny(1.0_dp) )             ) then
      call message()
      call message('***ERROR: L2_variable_init: Resolution of meteorology differs in basin: ', &
           trim(adjustl(num2str(iBasin))))
      stop
    end if

  
    ! cellfactor = leve1-2 / level-0
    cellFactor = level2%cellsize(iBasin) / level0%cellsize(iBasin)

    ! allocation and initalization of mask at level-2
    allocate( mask2(nrows2, ncols2) )
    mask2(:,:) = .FALSE.

    ! create mask at level-2
    do j = 1, ncols0
       jc = ceiling( real(j,dp)/cellFactor )
       do i = 1, nrows0
          if ( .NOT. mask0(i,j) ) cycle
          ic = ceiling( real(i,dp)/cellFactor )
          mask2(ic,jc) = .TRUE.
       end do
    end do
    
    ! no. of valid cells at level-2
    nCells2 = count( mask2 )

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    if (iBasin .eq. 1) then
 
       ! allocate space
       allocate(basin%L2_iStart     (nBasins))
       allocate(basin%L2_iEnd       (nBasins))
       allocate(basin%L2_iStartMask (nBasins))
       allocate(basin%L2_iEndMask   (nBasins))    

       ! basin information
       basin%L2_iStart(iBasin) = 1_i4
       basin%L2_iEnd  (iBasin) = basin%L2_iStart(iBasin) + nCells2 - 1_i4

       basin%L2_iStartMask(iBasin) = 1_i4
       basin%L2_iEndMask  (iBasin) = basin%L2_iStartMask(iBasin) + nrows2*ncols2 - 1_i4

    else

       ! basin information
       basin%L2_iStart(iBasin) = basin%L2_iEnd(iBasin-1) + 1_i4
       basin%L2_iEnd  (iBasin) = basin%L2_iStart(iBasin) + nCells2 - 1_i4

       basin%L2_iStartMask(iBasin) = basin%L2_iEndMask(iBasin-1) + 1_i4
       basin%L2_iEndMask  (iBasin) = basin%L2_iStartMask(iBasin) + nrows2*ncols2 - 1_i4

    end if

    call append( basin%L2_Mask,  RESHAPE( mask2, (/nrows2*ncols2/)  )  )

    ! free space
    deallocate(mask0, mask2)

  end subroutine L2_variable_init

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
  subroutine chunk_config( ii, tt, read_flag, start_date, end_date )
    !
    use mo_kind,             only: i4
    use mo_julian,           only: caldat
    use mo_global_variables, only: simPer
    use mo_mhm_constants,    only: nodata_dp
    !
    implicit none
    !
    ! input variables
    integer(i4), intent(in)  :: ii  ! basin id
    integer(i4), intent(in)  :: tt  ! current timestep
    !
    ! output variables
    logical,     intent(out) :: read_flag  ! indicate whether reading data should be read
    integer(i4), intent(out) :: start_date ! start date of read period
    integer(i4), intent(out) :: end_date   ! end date of read period
    !
    ! local variables
    integer(i4)              :: year     ! current year
    integer(i4)              :: month    ! current month
    integer(i4)              :: day      ! current day

    ! initialize
    read_flag  = .false.
    start_date = int( nodata_dp, i4 )
    end_date   = int( nodata_dp, i4 )

    ! evaluate date and timeStep_model_inputs to get read_flag -------
    read_flag = is_read( tt )
    !
    ! determine start and end date of chunk to read
    if ( read_flag ) call chunk_size( tt, start_date, end_date )
    print *, tt, read_flag, start_date, end_date
    !
  end subroutine chunk_config

  ! ------------------------------------------------------------------
  !
  ! function is_read
  !
  ! evaluates the current date and read condition
  !
  ! author: Stephan Thober
  !
  ! created: Jun 2014
  ! ------------------------------------------------------------------
  function is_read( tt )
    
    use mo_kind,             only: i4
    use mo_global_variables, only: simPer, timeStep_model_inputs, timestep
    use mo_message,          only: message
    use mo_julian,           only: caldat

    ! input variables
    integer(i4), intent(in) :: tt ! timestep
    
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
    Ndays        = ( tt * timestep ) / 24_i4 + 1
    Ndays_before = ( ( tt - 1_i4 ) * timestep ) / 24_i4 + 1
    
    ! evaluate cases of given timeStep_model_inputs
    select case( timeStep_model_inputs )
    case(0)  ! only at the beginning of the period
       if ( tt .eq. 1_i4 ) is_read = .true.
    case(1:) ! every timestep with frequency timeStep_model_inputs
       if ( mod( tt - 1_i4, 24 ) .eq. 0_i4 ) then
          if ( mod( (tt - 1_i4) / 24_i4 , timeStep_model_inputs ) .eq. 0_i4 ) is_read = .true.
       end if
    case(-1) ! every day
       if ( Ndays .ne. Ndays_before ) is_read = .true.
    case(-2) ! every month
       if ( Ndays .ne. Ndays_before ) then
          ! calculate months
          call caldat( simPer%julStart + Ndays, dd = day, mm = month, yy = year )
          call caldat( simPer%julStart + Ndays_before, dd = day_before, mm = month_before, yy = year_before )
          if ( month .ne. month_before ) is_read = .true.
       end if
    case(-3) ! every year
       if ( Ndays .ne. Ndays_before ) then
          ! calculate months
          call caldat( simPer%julStart + Ndays, dd = day, mm = month, yy = year )
          call caldat( simPer%julStart + Ndays_before, dd = day_before, mm = month_before, yy = year_before )
          if ( year .ne. year_before ) is_read = .true.
       end if
    case default ! not specified correctly
       call message('ERROR*** mo_meteo_forcings: function is_read: timStep_model_inputs not specified correctly!')
       stop
    end select
  
  end function is_read
  ! ------------------------------------------------------------------
  !
  ! subroutine chunk_size
  !
  ! calculating start and end date of chunk 
  !
  ! author: Stephan Thober
  !
  ! created: 2.6.2014
  ! 
  ! ------------------------------------------------------------------
  subroutine chunk_size( tt, start_date, end_date )
    
    use mo_kind,             only: i4
    use mo_global_variables, only: simPer, timeStep_model_inputs, timestep
    use mo_message,          only: message
    use mo_julian,           only: caldat, julday
    
    implicit none
    
    ! input variables
    integer(i4), intent(in)  :: tt
    
    ! output variables
    integer(i4), intent(out) :: start_date   ! julian start date
    integer(i4), intent(out) :: end_date     ! julian end date

    ! local variables
    integer(i4)              :: Ndays        ! number of simulated days
    integer(i4)              :: day          ! day
    integer(i4)              :: month        ! months
    integer(i4)              :: year         ! years

    ! calculate date of start date
    Ndays        = ( tt * timestep ) / 24_i4 + 1

    ! get start date
    start_date   = simPer%julStart + Ndays
    
    ! calculate length of the period
    select case ( timeStep_model_inputs )
    case(0)  ! length of chunk has to cover whole period
       end_date = simPer%julEnd 
    case(1:) ! every timestep with frequency timeStep_model_inputs
       end_date= start_date + timeStep_model_inputs - 1
    case(-1) ! every day
       end_date = start_date + 1_i4
    case(-2) ! every month
       ! calculate date
       call caldat( simPer%julStart + Ndays, dd = day, mm = month, yy = year )
       ! increment month
       if ( month .eq. 12 ) then
          month = 1
          year  = year + 1
       else
          month = month + 1
       end if
       end_date = julday( dd = 1, mm = month, yy = year ) - 1
    case(-3) ! every year
       ! calculate date
       call caldat( simPer%julStart + Ndays, dd = day, mm = month, yy = year )
       end_date = julday( dd = 31, mm = 12, yy = year )
    case default ! not specified correctly
       call message('ERROR*** mo_meteo_forcings: chunk_size: timStep_model_inputs not specified correctly!')
       stop
    end select   
    
    
  end subroutine chunk_size
  !
END MODULE mo_meteo_forcings
