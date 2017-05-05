!> \file mo_prepare_gridded_LAI.f90

!> \brief Prepare daily LAI fields (e.g., MODIS data) for mHM

!> \details Prepare daily LAI fields(e.g., MODIS data) for mHM

!> \authors John Craven & Rohini Kumar
!> \date Aug 2013

MODULE mo_prepare_gridded_LAI

  ! This module provides routines to read daily gridded LAI data.

  ! Written  John Craven & Rohini Kumar, August 2013
  ! Modified from mo_meteo_forcings

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: prepare_gridded_daily_LAI_data
  PUBLIC :: prepare_gridded_mean_monthly_LAI_data
 
  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !        prepare_gridded_daily_LAI_data
  
  !     PURPOSE
  !>        \brief Prepare gridded daily LAI data 

  !>        \details Prepare gridded daily LAI data at Level-0 (e.g., using MODIS datasets)

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
  !>        \author John Craven & Rohini Kumar
  !>        \date Aug 2013
  !               Modified Matthias Cuntz & Juliane Mai, Nov 2014 - use meteo reading routines
  !
  subroutine prepare_gridded_daily_LAI_data(iBasin)
    
    use mo_global_variables,           only: dirgridded_LAI, inputFormat_gridded_LAI, &      
                                             simPer, L0_gridded_LAI, timeStep_LAI_input, L1_gridded_LAI
    use mo_init_states,                only: get_basin_info            ! get basin information
    use mo_append,                     only: append                    ! append vector
    use mo_read_meteo,                 only: read_meteo_bin            ! Read binary files
    use mo_read_forcing_nc,            only: read_forcing_nc           ! Read netCDF files
    
    use mo_global_variables,           only: readPer, level1, level2
    use mo_mhm_constants,              only: nodata_dp
    !use mo_init_states,                only: get_basin_info
    !use mo_read_meteo,                 only: read_meteo_bin
    !use mo_read_forcing_nc,            only: read_forcing_nc
    use mo_spatial_agg_disagg_forcing, only: spatial_aggregation, spatial_disaggregation
    !use mo_append,                     only: append                    ! append vector
    
    implicit none

    integer(i4),                          intent(in)    :: iBasin        ! Basin Id
    !character(len=*),                     intent(in)    :: dataPath      ! Data path
    !character(len=*),                     intent(in)    :: inputFormat   ! either 'bin' or 'nc'
    !real(dp), dimension(:,:),allocatable, intent(inout) :: dataOut1      ! Packed meteorological variable
    !real(dp),                   optional, intent(in)    :: lower         ! lower bound for data points
    !real(dp),                   optional, intent(in)    :: upper         ! upper bound for data points
    !character(len=*),           optional, intent(in)    :: ncvarName     ! name of the variable (for .nc files)

    
    integer(i4)                                :: nrows1, ncols1
    logical, dimension(:,:), allocatable       :: mask1
    integer(i4)                                :: ncells1

    integer(i4)                                :: nrows2, ncols2
    logical, dimension(:,:), allocatable       :: mask2

    !real(dp), dimension(:,:,:), allocatable    :: L2_data            ! meteo data at level-2 
    !real(dp), dimension(:,:,:), allocatable    :: L1_data            ! meteo data at level-1
    !real(dp), dimension(:,:), allocatable      :: L1_data_packed     ! packed meteo data at level-1 from 3D to 2D

    integer(i4)                                :: nTimeSteps
    real(dp)                                   :: cellFactorHbyM   ! level-1_resolution/level-2_resolution
    integer(i4)                                :: t
    integer(i4)                                :: k

                                     
    ! input 
    !integer(i4),                   intent(in)  :: iBasin    ! Basin Id
    !local variables
    ! basin info.
    integer(i4)                                :: nrows0, ncols0
    logical, dimension(:,:), allocatable       :: mask0
    integer(i4)                                :: ncells0
    !
    real(dp), dimension(:,:,:), allocatable    :: LAI1_3D     !data at level-1 [nRow X nCols X nTimeSteps]
    real(dp), dimension(:,:), allocatable    :: LAI1_2D     !data at level-1 [nRow X nCols X nTimeSteps]

    real(dp), dimension(:,:,:), allocatable    :: LAI0_3D     !data at level-0 [nRow X nCols X nTimeSteps]
    real(dp), dimension(:,:), allocatable      :: LAI0_2D     !data at level-0 [nCells X nTimeSteps]

    !integer(i4)                                :: nTimeSteps
    !integer(i4)                                :: t

    ! get basic basin information at level-0
    call get_basin_info( iBasin, 0, nRows0, nCols0, nCells=nCells0, mask=mask0 )
    
    ! get basic basin information at level-1
    call get_basin_info( iBasin, 1, nrows1, ncols1, nCells=nCells1, mask=mask1 )
    
    ! make  basic basin information at level-2
    call get_basin_info( iBasin, 2, nrows2, ncols2, mask=mask2 )
      

    ! select case depending on input data format
    SELECT CASE( trim(inputFormat_gridded_LAI) )       

       ! netcdf file input option
       CASE('nc')
          CALL read_forcing_nc( dirgridded_LAI(iBasin), nRows0, nCols0, simPer(iBasin), &
               'lai', LAI0_3D, mask0, lower=0.0_dp, upper=30.0_dp, nctimestep=timeStep_LAI_input)
        !CUNEYD: Lai.nc is read at L0 level. We can upscale it to L1 and use it in Kc correction of PET  
        !CALL read_forcing_nc( dataPath, nRows2, nCols2, readPer, ncvarName, L2_data, mask2, &
        !       lower=lower, upper=upper )
        call spatial_aggregation(LAI0_3D, level2%cellsize(iBasin), level1%cellsize(iBasin), mask1, mask2, LAI1_3D)
       
       ! bin file input option
       CASE('bin')
          CALL read_meteo_bin( dirgridded_LAI(iBasin), nRows0, nCols0, simPer(iBasin), &
               LAI0_3D, mask0, lower=0.0_dp, upper=30.0_dp)
       CASE DEFAULT
           stop '***ERROR: Not recognized input format'
           
    END SELECT

    ! pack variables
    nTimeSteps = size(LAI0_3D, 3)
    allocate( LAI0_2D(nCells0, nTimeSteps))

    do t = 1, nTimeSteps
       LAI0_2D(:,t) = pack( LAI0_3D(:,:,t), MASK=mask0(:,:) ) 
    end do
    
    ! append to Global variable
    call append(L0_gridded_LAI, LAI0_2D(:,:) )
    
    !free space
    deallocate(LAI0_2D, LAI0_3D) 
    
    
    
    ! pack variables
    nTimeSteps = size(LAI1_3D, 3)
    allocate( LAI1_2D(nCells1, nTimeSteps))

    do t = 1, nTimeSteps
       LAI1_2D(:,t) = pack( LAI1_3D(:,:,t), MASK=mask1(:,:) ) 
    end do
    
    print*,"LAI1 ", LAI1_2D(1,1)
    ! append to Global variable
    call append(L1_gridded_LAI, LAI1_2D(:,:) )
    
    !free space
    deallocate(LAI1_2D, LAI1_3D) 
 
 
end subroutine prepare_gridded_daily_LAI_data

  ! ------------------------------------------------------------------

  !     NAME
  !        prepare_gridded_mean_monthly_LAI_data
  
  !     PURPOSE
  !>        \brief prepare_gridded_mean_monthly_LAI_data

  !>        \details Prepare long term mean monthly gridded LAI data at Level-0 (e.g., using MODIS datasets)

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
  !>        \date Dec 2016
  !
  subroutine prepare_gridded_mean_monthly_LAI_data(iBasin)
    
    use mo_global_variables,           only: dirgridded_LAI, L0_gridded_LAI
    use mo_init_states,                only: get_basin_info            ! get basin information
    use mo_append,                     only: append                    ! append vector
    use mo_ncread,                     only: Get_NcDim, Get_NcVar, Get_NcVarAtt
    use mo_message,                    only: message
    use mo_string_utils,               only: num2str
    use mo_utils,                      only: eq
                                           
    implicit none
    ! input 
    integer(i4),                   intent(in)  :: iBasin    ! Basin Id
    !local variables
    ! basin info.
    integer(i4)                                :: nrows0, ncols0
    logical, dimension(:,:), allocatable       :: mask0
    integer(i4)                                :: ncells0
    !
    real(dp), dimension(:,:,:), allocatable    :: LAI0_3D     !data at level-0 [nRow X nCols X nTimeSteps]
    real(dp), dimension(:,:), allocatable      :: LAI0_2D     !data at level-0 [nCells X nTimeSteps]

    integer(i4)                                :: nTimeSteps
    integer(i4)                                :: t

    !
    character(256)                             :: fName        ! name of NetCDF file
    character(256)                             :: AttValues    ! netcdf attribute values
    integer(i4)                                :: datatype     ! datatype of attribute
    integer(i4), dimension(5)                  :: dimen        ! dimension for NetCDF file
    real(dp)                                   :: nodata_value ! data nodata value
    
    
    
    ! get basic basin information at level-0
    call get_basin_info( iBasin, 0, nRows0, nCols0, nCells=nCells0, mask=mask0 )

    fName = trim( dirgridded_LAI(iBasin) ) // trim('lai.nc')
    
    ! get dimensions
    dimen = Get_NcDim( trim(fName), 'lai' )
    if ( (dimen(1) .ne. nRows0) .or. (dimen(2) .ne. nCols0) ) then
       stop '***ERROR: read_forcing_nc: mHM generated x and y are not matching NetCDF dimensions'
    end if
    if ( dimen(3) .ne. 12 ) then
       stop '***ERROR: read_forcing_nc: the time dimension of LAI NetCDF file under the option-1 is not 12'
    end if


    ! determine no data value
    call Get_NcVarAtt( trim(fName), 'lai', '_FillValue', AttValues, dtype=datatype)
    ! convert to number
    read(AttValues, *) nodata_value
        
    call Get_NcVar( trim(fName), 'lai', LAI0_3D )

    ! start checking values
    do t = 1, dimen(3) 
       ! checking for naodata values if optional nocheck is given
       if (any(eq(LAI0_3D(:,:,t),nodata_value) .and. (mask0))) then
          call message('***ERROR: read_forcing_nc: nodata value within basin ')
          call message('          boundary in variable: ', 'lai')
          call message('          at timestep         : ', trim(num2str(t)))
          stop
       end if
       ! optional check
       if ( any( (LAI0_3D(:,:,t) .lt. 0.0_dp) .AND. mask0(:,:) )  ) then
          call message('***ERROR: read_forcing_nc: values in variable lai are lower than ', trim(num2str(0,'(F7.2)')) )
          call message('          at timestep  : ', trim(num2str(t)))
          call message('File: ', trim(fName))
          call message('Minval at timestep: ', trim(num2str(minval(LAI0_3D(:,:,t)),'(F7.2)')))
          call message('Total minval: ', trim(num2str(minval(LAI0_3D(:,:,:)),'(F7.2)')))
          stop
       end if

       if ( any( (LAI0_3D(:,:,t) .gt. 30.0_dp) .AND. mask0(:,:) )  ) then
          call message('***ERROR: read_forcing_nc: values in variable lai are greater than ', trim(num2str(30,'(F7.2)')) )
          call message('          at timestep  : ', trim(num2str(t)))
          call message('File: ', trim(fName))
          call message('Maxval at timestep: ', trim(num2str(maxval(LAI0_3D(:,:,t)),'(F7.2)')))
          call message('Total maxval: ', trim(num2str(maxval(LAI0_3D(:,:,:)),'(F7.2)')))
          stop
       end if
     end do
    
    ! pack variables
    nTimeSteps = size(LAI0_3D, 3)
    allocate( LAI0_2D(nCells0, nTimeSteps))
    do t = 1, nTimeSteps
       LAI0_2D(:,t) = pack( LAI0_3D(:,:,t), MASK=mask0(:,:) ) 
    end do
    
    ! append to Global variable
    call append(L0_gridded_LAI, LAI0_2D(:,:) )
    
    !free space
    deallocate(LAI0_2D, LAI0_3D) 
 
 
end subroutine prepare_gridded_mean_monthly_LAI_data


END MODULE mo_prepare_gridded_LAI
