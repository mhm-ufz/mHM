!> \file mo_spatial_agg_disagg_forcing.f90

!> \brief Spatial aggegation or disaggregation of meteorological input data.

!> \details This module contains two subroutines to upscale and downscale, respectively,\n
!>          the level-2 meterological inputs to a required Level-1 hydrological spatial \n
!>          resolution. 

!> \authors Rohini Kumar
!> \date Jan 2013

MODULE mo_spatial_agg_disagg_forcing

  ! This module provides routines for spatial aggegation or disaggregation of meteorological input data.

  ! Written  Rohini Kumar, Jan 2013
  ! Modified 

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: spatial_aggregation      ! Spatial aggregation or upscaling
  PUBLIC :: spatial_disaggregation   ! Spatial disaggregation or downscaling

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         spatial_aggregation

  !     PURPOSE
  !>        \brief Spatial aggregation of meterological variables

  !>        \details Aggregate (or upscale) the given level-2 meteorological data to the \n 
  !>                 required level-1 spatial resolution for the mHM run. 

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:,:,:)   :: data2"         Level-2 meteorological data \n
  !>                                                                   dim1=nRows, dim2=nCols, dim3=nTimeSteps
  !>        \param[in] "real(dp)                     :: cellsize2"     Level-2 resolution 
  !>        \param[in] "real(dp)                     :: cellsize1"     Level-1 resolution 
  !>        \param[in] "logical, dimension(:,:)      :: mask1"         Level-1 basin mask \n
  !>        \param[in] "logical, dimension(:,:)      :: mask2"         Level-2 basin mask \n
  !>                                                                   dim1=nRows, dim2=nCols

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), allocatable, dimension(:,:,:)  :: data1" 
  !>           return level-1 meteorological data; DIMENSION [nRows_L1, nCols_L1, nTimeSteps] 

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
  !          Modified Rohini Kumar,   Nov 2013 - data1 changed from intent(inout) to intent(out)
  !          Modified RK, MZ, DS,     May 2014 - added mask2

  subroutine spatial_aggregation(data2, cellsize2, cellsize1, mask1, mask2, data1)
    
    use mo_mhm_constants, only: nodata_dp

    implicit none

    real(dp), dimension(:,:,:),              intent(in)  :: data2       ! Level-2 data
    real(dp),                                intent(in)  :: cellsize2   ! Level-2 resolution  
    real(dp),                                intent(in)  :: cellsize1   ! Level-1 resolution 
    logical, dimension(:,:),                 intent(in)  :: mask1       ! Level-1 mask
    logical, dimension(:,:),                 intent(in)  :: mask2       ! Level-2 mask

    real(dp), dimension(:,:,:), allocatable, intent(out) :: data1       ! Level-1 data
        
    ! local variables
    integer(i4)                                :: nr2, nc2  ! No. of rows and cols at Level-1 
    integer(i4)                                :: nr1, nc1  ! No. of rows and cols at Level-1 
    real(dp)                                   :: cellFactor 
    integer(i4), dimension(:,:), allocatable   :: nTCells
    integer(i4)                                :: nTimeSteps
    integer(i4)                                :: i, j, ic, jc, t

    ! get number of rows and cols at level-2 from mask2
    ! and the total time steps 
    nr2        = size( data2, 1 )
    nc2        = size( data2, 2 )
    nTimeSteps = size( data2, 3 )
       
    ! get number of rows and cols at level-1 from mask1
    nr1 = size( mask1, 1 )
    nc1 = size( mask1, 2 )

    !-----------------------------------------------------------------------
    ! Allocate and initalize nTCells which comprises
    ! of number of L2 cells that belongs to a given L1 cell
    ! NOTE:: 1) cell size of L1 > L2 (see CellFactor)
    !        2) nTCells is estimated over the valid masked domain only
    !-----------------------------------------------------------------------
    
    ! cellFactor = level-1 resolution (hydro) / level-2 resolution (meteo)
    cellFactor = cellsize1 / cellsize2
    
    ! nTCells calculations
    allocate( nTCells(nr1, nc1) )
    nTCells(:,:) = 0

    do j=1, nc2
       jc = ceiling( real(j,dp)/cellFactor )
       do i=1, nr2
          ic = ceiling(real(i,dp)/cellFactor)
          if( .not. mask2(i,j) ) cycle
          nTCells(ic,jc) = nTcells(ic,jc) + 1  
       end do
    end do
    
         
    ! allocate and initalize L1_data
    allocate( data1(nr1, nc1, nTimeSteps) )
    data1(:,:,:) = 0.0_dp
   
    ! time loop
    do t = 1, nTimeSteps
    
       ! perform spatial aggregation
       do j=1, nc2
          jc = ceiling( real(j,dp)/cellFactor )
          do i=1, nr2
             ic = ceiling(real(i,dp)/cellFactor)
             
             ! only in valid masked area
             if( .not. mask2(i,j) ) cycle
             data1(ic,jc,t) = data1(ic,jc,t) + data2(i,j, t)
             
          end do
       end do
       
       ! perform spatial average only over valid masked domain
       ! out of the masked domain nTCells(:,:) = 0
       where( mask1 )
          data1(:,:,t) = data1(:,:,t) / real( nTcells(:,:), dp )
       elsewhere
          data1(:,:,t) = nodata_dp
       endwhere

    end do
    
    ! free space 
    deallocate(nTCells)

  end subroutine spatial_aggregation
  
  
  ! ------------------------------------------------------------------

  !     NAME
  !         spatial_disaggregation

  !     PURPOSE
  !>        \brief Spatial disaggregation of meterological variables

  !>        \details Disaggregate (or downscale) the given level-2 meteorological data to the \n 
  !>                 required level-1 spatial resolution for the mHM run. 

  !     CALLING SEQUENCE

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:,:,:)   :: data2"         Level-2 meteorological data \n
  !>                                                                   dim1=nRows, dim2=nCols, dim3=nTimeSteps
  !>        \param[in] "real(dp)                     :: cellsize2"     Level-2 resolution 
  !>        \param[in] "real(dp)                     :: cellsize1"     Level-1 resolution 
  !>        \param[in] "logical, dimension(:,:)      :: mask1"         Level-1 basin mask \n
  !>        \param[in] "logical, dimension(:,:)      :: mask2"         Level-2 basin mask \n
  !>                                                                   dim1=nRows, dim2=nCols

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), allocatable, dimension(:,:,:)  :: data1" 
  !>           return level-1 meteorological data; DIMENSION [nRows_L1, nCols_L1, nTimeSteps] 

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
  !          Modified Rohini Kumar,   Nov 2013 - data1 changed from intent(inout) to intent(out)
  !          Modified RK, MZ, DS,     May 2014 - added mask2

  subroutine spatial_disaggregation(data2, cellsize2, cellsize1, mask1, mask2, data1)

    use mo_mhm_constants, only: nodata_dp

    implicit none

    real(dp), dimension(:,:,:),                 intent(in) :: data2       ! Level-2 data
    real(dp),                                   intent(in) :: cellsize2   ! Level-2 resolution  
    real(dp),                                   intent(in) :: cellsize1   ! Level-1 resolution 
    logical, dimension(:,:),                    intent(in) :: mask1       ! Level-1 mask
    logical, dimension(:,:),                    intent(in) :: mask2       ! Level-2 mask

    real(dp), dimension(:,:,:), allocatable, intent(out)   :: data1       ! Level-1 data
        
    ! local variables
    integer(i4)                                :: nr1, nc1                ! No. of rows and cols at Level-1 

    real(dp)                                   :: cellFactor 
    integer(i4)                                :: nTimeSteps
    integer(i4)                                :: i, j, t, ic, jc

    ! get number of rows and cols at level-2 from mask2
    nr1 = size( mask1, 1)
    nc1 = size( mask1, 2)

    ! cellFactor = level-2 resolution (meteo) / level-1 resolution (hydro)
    cellFactor = cellsize2 / cellsize1
    
    ! total time steps 
    nTimeSteps = size(data2, 3)

    ! allocate and initalize L1_data
    allocate( data1(nr1, nc1, nTimeSteps) )
    data1(:,:,:) = nodata_dp
    
    ! over the time loop
    do t = 1, nTimeSteps
    
      ! spatial disaggregation
      do j=1, nc1
         jc = ceiling( real(j,dp)/cellFactor )
         do i=1, nr1
           ic = ceiling(real(i,dp)/cellFactor)
           ! only over the valid masked area
           if( .not. mask2(ic,jc) ) cycle
           data1(i,j,t) = data2(ic,jc, t)
          end do
      end do
      
    end do    
   
  end subroutine spatial_disaggregation


END MODULE mo_spatial_agg_disagg_forcing
