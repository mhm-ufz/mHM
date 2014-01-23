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
  subroutine prepare_meteo_forcings_data(iBasin)
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
    cellFactorHbyM = level1%cellsize / level2%cellsize 

    ! upscaling & packing
    if(cellFactorHbyM .gt. 1.0_dp) then 
        call spatial_aggregation(L2_data, level2%cellsize, level1%cellsize, mask1, L1_data)
    ! downscaling   
    elseif(cellFactorHbyM .lt. 1.0_dp) then
        call spatial_disaggregation(L2_data, level2%cellsize, level1%cellsize, mask1, L1_data)
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
       allocate( level2%nrows     (nBasins) )
       allocate( level2%ncols     (nBasins) )
       allocate( level2%xllcorner (nBasins) )
       allocate( level2%yllcorner (nBasins) )
     end if

    ! read header file 
    ! NOTE: assuming the header file for all metero variables are same as that of precip.
    !       A counter check for this assumption is perfromed in the read_meteo_bin file 
    
    fName =  trim(adjustl(dirPrecipitation(iBasin))) // trim(adjustl(file_meteo_header))
    call read_header_ascii( trim(fName), umeteo_header,   &
                            level2%nrows(iBasin), level2%ncols(iBasin), level2%xllcorner(iBasin), &
                            level2%yllcorner(iBasin), level2%cellsize,  level2%nodata_value       )
  
   ! level-0 information
   call get_basin_info( iBasin, 0, nrows0, ncols0, mask=mask0,                         &
                        xllcorner=xllcorner0, yllcorner=yllcorner0, cellsize=cellsize0 ) 
   ! grid information
   call calculate_grid_properties( nrows0, ncols0, xllcorner0, yllcorner0, cellsize0, nodata_dp,          &
                                   level2%cellsize, &
                                   nrows2, ncols2, xllcorner2, yllcorner2, cellsize2,level2%nodata_value )
   ! check
   if (  (ncols2     .ne.  level2%ncols(iBasin))         .or. &
         (nrows2     .ne.  level2%nrows(iBasin))         .or. &
         ( abs(xllcorner2 - level2%xllcorner(iBasin)) .gt. tiny(1.0_dp) )     .or. &
         ( abs(yllcorner2 - level2%yllcorner(iBasin)) .gt. tiny(1.0_dp) )     .or. &
         ( abs(cellsize2  - level2%cellsize)         .gt. tiny(1.0_dp) )             ) then
      call message()
      call message('***ERROR: L2_variable_init: Resolution of meteorology differs in basin: ', &
           trim(adjustl(num2str(iBasin))))
      stop
    end if

  
    ! cellfactor = leve1-2 / level-0
    cellFactor = level2%cellsize / level0%cellsize

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

END MODULE mo_meteo_forcings
