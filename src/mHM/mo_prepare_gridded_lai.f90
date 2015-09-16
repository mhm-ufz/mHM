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
                                             simPer, L0_gridded_LAI, timeStep_LAI_input
    use mo_init_states,                only: get_basin_info            ! get basin information
    use mo_append,                     only: append                    ! append vector
    use mo_read_meteo,                 only: read_meteo_bin            ! Read binary files
    use mo_read_forcing_nc,            only: read_forcing_nc           ! Read netCDF files
                                           
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

    ! get basic basin information at level-0
    call get_basin_info( iBasin, 0, nRows0, nCols0, nCells=nCells0, mask=mask0 )

    ! select case depending on input data format
    SELECT CASE( trim(inputFormat_gridded_LAI) )       

       ! netcdf file input option
       CASE('nc')
          CALL read_forcing_nc( dirgridded_LAI(iBasin), nRows0, nCols0, simPer(iBasin), &
               'lai', LAI0_3D, mask0, lower=0.0_dp, upper=30.0_dp, nctimestep=timeStep_LAI_input)
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
 
 
end subroutine prepare_gridded_daily_LAI_data



END MODULE mo_prepare_gridded_LAI
