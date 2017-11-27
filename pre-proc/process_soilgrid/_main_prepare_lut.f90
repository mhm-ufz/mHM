!> \file main.f90
!**********************************************************************************
!> PREPERATION OF THE LUT \n
!>  PURPOSE: \n  
!
!>  AUTHOR: \n   
!>            Rohini Kumar, UFZ, March, 2016 \n
!
!  HISTORY
!         Created   R. Kumar  March, 2016          main structure
!
!**********************************************************************************
 Program main
  use mo_kind
  implicit none
  !
  ! soil related properties
  integer(i4), parameter               :: nsoil_layers  = 6
  integer(i4), dimension(nsoil_layers) :: soil_up_depth  = (/0,  50,  150, 300,  600, 1000/)
  integer(i4), dimension(nsoil_layers) :: soil_low_depth = (/50, 150, 300, 600, 1000, 2000/)
  !
  ! grid header info.
  integer(i4)                                :: ncols               ! number of columns
  integer(i4)                                :: nrows               ! number of rows
  real(sp)                                   :: xllcorner           ! x coordinate of the lowerleft corner
  real(sp)                                   :: yllcorner           ! y coordinate of the lowerleft corner
  real(sp)                                   :: cellsize            ! cellsize x = cellsize y
  integer(i4)                                :: nodata_value        ! code to define the mask
  integer(i4), dimension(:,:,:), allocatable :: sd, cl, bd          ! [nrows, ncols, nlayers]
  ! cell level
  integer(i4)                                :: ncells              ! total number of cells = no. of soil types
  integer(i4), dimension(:,:), allocatable   :: soilid
  !
  character(256)                             :: dummy, dataPath_in, dataPath_out, fName
  integer(i4)                                :: ii, jj, kk, ll
  !
  ! datapath .....
  dataPath_in = "/home/rkumar/projects/edge/soil_process/ascii/"
  dataPath_out = "/data/edge/data/processed/morph/process/"
  !
  ! read gridded files
  do kk = 1, nsoil_layers
     print*, 'READING LAYER', kk
     !
     write(dummy,100) kk
100  format(i2.2,'.txt')
     !
     fName = trim(dataPath_in)//trim('bd')//trim(dummy) 
     open(10, file=trim(fName), status='old', action = 'read')
     fName = trim(dataPath_in)//trim('cl')//trim(dummy) 
     open(20, file=trim(fName), status='old', action = 'read')
     fName = trim(dataPath_in)//trim('sn')//trim(dummy) 
     open(30, file=trim(fName), status='old', action = 'read')
     !
     ! only once -- all files must have same header
     if( kk == 1 ) then
        read(10,*) dummy, ncols
        read(10,*) dummy, nrows
        read(10,*) dummy, xllcorner
        read(10,*) dummy, yllcorner
        read(10,*) dummy, cellsize
        read(10,*) dummy, nodata_value
        allocate( bd(nrows, ncols, nsoil_layers) )
        allocate( cl(nrows, ncols, nsoil_layers) )
        allocate( sd(nrows, ncols, nsoil_layers) )
        bd(:,:,:) = -9999
        cl(:,:,:) = -9999
        sd(:,:,:) = -9999
        ! to the start of the file
        rewind(10)
     end if
     ! read the headers and the gridded data
     do ii = 1, 6
        read(10,*) dummy
        read(20,*) dummy
        read(30,*) dummy
     end do
     !
     do ii = 1, nrows
        read(10,*) ( bd(ii,jj,kk), jj=1,ncols )
        read(20,*) ( cl(ii,jj,kk), jj=1,ncols )
        read(30,*) ( sd(ii,jj,kk), jj=1,ncols )
     end do
     close(10); close(20); close(30)
     !
  end do
  !
  !-------------------------------------------------------------------------------------------
  ! READ AND PROCESS MASKED DATA SETS
  !         in mHM FORMAT
  !-------------------------------------------------------------------------------------------
  !!
  !! check that all three soil properties are present in all layers
  !! i.e. all have valid data
  print*, 'nvalid_bd', 'nbd_clay', 'nbd_sand'
  print*, count( bd(:,:,:) .ge. 0 ), count( cl(:,:,:) .ge. 0 ), count( sd(:,:,:) .ge. 0 )
  !!
  !! remove zero values
  where( bd(:,:,:) .eq. 0 ) bd(:,:,:) = 1500
  where( cl(:,:,:) .eq. 0 ) cl(:,:,:) = 33
  where( sd(:,:,:) .eq. 0 ) sd(:,:,:) = 33
  !!
  !! potential number of soil types
  ncells = count( sd(:,:,1) .ge. 0 )
  allocate( soilid(nrows, ncols) )
  soilid(:,:) = -9999
  !!
  !! write results -- LUT
  fName = trim(dataPath_out)//trim('soil_lut.txt') 
  open(40, file = trim(fName), status = 'unknown', action = 'write')
  write(40,400)'nSoil_Types', (ncells+1)
  write(40,401)'ID', 'HORIZON', 'UD[mm]','LD[mm]','CLAY[%]','SAND[%]', 'Bd[gcm-3]'
  !
  kk = 0
  do jj = 1, ncols
     do ii = 1, nrows
        if( sd(ii,jj,1) .lt. 0 ) cycle
        kk = kk + 1
        soilid(ii,jj) = kk
        do ll = 1, nsoil_layers
           write(40, 402) soilid(ii,jj), ll, soil_up_depth(ll),  soil_low_depth(ll), &
                          cl(ii,jj,ll), sd(ii,jj,ll), real(bd(ii,jj,ll),sp)/1000.0_sp
        end do
     end do
  end do
  !!
  !! add one extra soil type to fill missing estimates
  where( soilid(:,:) < 0 ) soilid(:,:) = ncells+1
  write(40, 402) (ncells+1), 1,   0,  300, 33, 33, 1.5_sp
  write(40, 402) (ncells+1), 2, 300, 2000, 33, 33, 1.5_sp
  close(40)
!!!
!!! formats
  400 format(A15, I15)
  401 format(A15, 6A15)
  402 format(I10, 1X, I3, 2I5, 2I4, f6.3)
  !!
  !! WRITE soil id file
  fName = trim(dataPath_out)//trim('soil_id.asc') 
  open (50, file= trim(fName), status='unknown', action = 'write')
  write(50,1)'ncols       ', ncols
  write(50,1)'nrows       ', nrows
  write(50,2)'xllcorner   ', xllcorner
  write(50,2)'yllcorner   ', yllcorner
  write(50,2)'cellsize    ', cellsize
  write(50,1)'NODATA_value', nodata_value
  do ii=1,nrows
     write(50, 981) ( soilid(ii,jj), jj=1,ncols )
  end do
  close(50)
!!!
!!! formats
  1 format( a12, 2x, i15   )
  2 format( a12, 2x, f15.6 )
981 format( 1000000(i10, 1x) )




  !
 end program main




