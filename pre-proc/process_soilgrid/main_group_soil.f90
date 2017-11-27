!> \file main.f90
!**********************************************************************************
!> PROCESSOR OF SOIL GRID DTASETS \n
!>  PURPOSE: \n  
!>           Convert the soil grid datasest to mHM format
!>           Get gridded soil textural dataset from www.soilgrids.org
!>           There should be six horizon informations from the soil grid
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
   use mo_moment, only: average, stddev
  implicit none
  !!
  !! soil related properties
  integer(i4)                                :: nsoil_types, nsoil_types_reduced              ! total no. of soil types
  integer(i4), parameter                     :: nsoil_layers   = 6
  integer(i4), dimension(nsoil_layers)       :: soil_up_depth  = (/0,  50,  150, 300,  600, 1000/)
  integer(i4), dimension(nsoil_layers)       :: soil_low_depth = (/50, 150, 300, 600, 1000, 2000/)
  !!
  !! grid header info.
  integer(i4)                                :: ncols               ! number of columns
  integer(i4)                                :: nrows               ! number of rows
  real(sp)                                   :: xllcorner           ! x coordinate of the lowerleft corner
  real(sp)                                   :: yllcorner           ! y coordinate of the lowerleft corner
  real(sp)                                   :: cellsize            ! cellsize x = cellsize y
  integer(i4)                                :: nodata_value        ! code to define the mask
  type soil_prop_original
     integer(i4), dimension(:), allocatable  :: sand, clay
     real(sp), dimension(:), allocatable     :: bd
  end type soil_prop_original
  type(soil_prop_original), dimension(:), allocatable :: soil_orig
  integer(i4), dimension(:,:), allocatable   :: soilid, soilid_horizon
  !!
  !! output
  integer(i4)                                :: ncols_out               ! number of columns
  integer(i4)                                :: nrows_out               ! number of rows
  real(sp)                                   :: xllcorner_out           ! x coordinate of the lowerleft corner
  real(sp)                                   :: yllcorner_out           ! y coordinate of the lowerleft corner
  real(sp)                                   :: cellsize_out            ! cellsize x = cellsize y
  integer(i4)                                :: nodata_value_out        ! code to define the mask
  type soil_prop_reduced
     integer(i4)                             :: soil_id
     integer(i4)                             :: nsoils                  ! frequency 
     real(sp)                                :: bd_mean, bd_stdev       ! bulk density
     real(sp), dimension(:), allocatable     :: bd                      ! consolidate all bd for this soil grid
  end type soil_prop_reduced
  type(soil_prop_reduced), dimension(:,:), allocatable :: grid_out
  !!
  character(256)                             :: dummy, dataPath, fName
  integer(i4)                                :: ii, jj, kk, ll, id, hh
  integer(i4)                                :: sand, clay
  real(sp)                                   :: bulk_density

  !!
  !! program starts here
  dataPath = "/data/edge/data/processed/mhm_input/morph/"

  !!
  !! assign values (to plot) and allocate variables
  nrows_out = 100   !>>> sand content in 100 parts
  ncols_out = 100   !>>> clay content in 100 parts
  xllcorner_out = 0.0_sp
  yllcorner_out = 100.0_sp
  cellsize_out  = 1.0_sp
  nodata_value_out = -9999.0_sp
  allocate( grid_out(nrows_out, ncols_out) )
  grid_out(:,:)%nsoils = 0
  !!
  !! read the original soil LUT to get the frequency grid
  fName = trim(dataPath)//trim('soil_classdefinition.txt') 
  open(10, file = trim(fName), status = 'old', action = 'read')
  read(10,*) dummy, nsoil_types
  read(10,*) dummy
  allocate( soil_orig(nsoil_types) )
  do kk = 1, (nsoil_types-1)*nsoil_layers  !! last soil is default and had only two layers
     read(10,*) id, hh, dummy, dummy, clay, sand, bulk_density
     !! limits
     if( clay == 0 ) clay = 1
     if( sand == 0 ) sand = 1
     if (hh == 1) then
        allocate( soil_orig(id)%sand(nsoil_layers) )
        allocate( soil_orig(id)%clay(nsoil_layers) )
        allocate( soil_orig(id)%bd  (nsoil_layers) )
     end if
     soil_orig(id)%sand(hh) = sand
     soil_orig(id)%clay(hh) = clay     
     soil_orig(id)%bd  (hh) = bulk_density
     grid_out(sand,clay)%nsoils = grid_out(sand,clay)%nsoils + 1
  end do
  close(10)
  !! last soil layer info....
  id = nsoil_types
  allocate( soil_orig(id)%sand(nsoil_layers) )
  allocate( soil_orig(id)%clay(nsoil_layers) )
  allocate( soil_orig(id)%bd  (nsoil_layers) )
  soil_orig(id)%sand(:) = 33
  soil_orig(id)%clay(:) = 33     
  soil_orig(id)%bd  (:) = 1.5_sp 
  grid_out(33,33)%nsoils = grid_out(33,33)%nsoils + nsoil_layers
  !!
  print*, 'reading of the soil database over'
  !!
  !! make unique reduced soil class ids
  nsoil_types_reduced = 0
  grid_out%soil_id    = int(nodata_value_out)
  do jj = 1, ncols_out
     do ii = 1, nrows_out
        if( grid_out(ii,jj)%nsoils < 1 ) cycle
         nsoil_types_reduced     = nsoil_types_reduced + 1
         grid_out(ii,jj)%soil_id = nsoil_types_reduced
     end do
  end do
  print*, 'No of reduced soil class', nsoil_types_reduced

  !!===================================================================================
  !!
  !! read soil id
  fName = trim(dataPath)//trim('soil_class.asc') 
  open(10, file=trim(fName), status='old', action='read')
  read(10,*) dummy, ncols
  read(10,*) dummy, nrows
  read(10,*) dummy, xllcorner
  read(10,*) dummy, yllcorner
  read(10,*) dummy, cellsize
  read(10,*) dummy, nodata_value
  allocate( soilid(nrows, ncols) )
  do ii = 1, nrows
     read(10,*) ( soilid(ii,jj), jj=1,ncols )
  end do
  close(10)
  print*, 'reading of the original soil-id over'
  !!
  !! process horizon wise new information
  allocate( soilid_horizon(nrows, ncols) )
  do hh = 1, nsoil_layers
     soilid_horizon(:,:) = int(nodata_value)
     do jj = 1, ncols
        do ii = 1, nrows
           id = soilid(ii,jj) 
           if( id < 0 ) cycle
           sand = soil_orig(id)%sand(hh)
           clay = soil_orig(id)%clay(hh)          
           soilid_horizon(ii,jj) = grid_out(sand,clay)%soil_id
        end do
     end do
     !!     
     !! write asc files
     write(fName,100) soil_up_depth(hh), soil_low_depth(hh)
100  format( 'soil_class_horizon_depth_', i5.5, '-',i5.5, 'mm.asc')
     fName = trim(dataPath)//trim(fName) 
     open(50, file=trim(fName), status='unknown', action = 'write')
     write(50,1) 'ncols       ', ncols
     write(50,1) 'nrows       ', nrows
     write(50,2) 'xllcorner   ', xllcorner
     write(50,2) 'yllcorner   ', yllcorner
     write(50,2) 'cellsize    ', cellsize
     write(50,1) 'NODATA_value', int(nodata_value)
     do ii=1,nrows
        write(50, 980) ( soilid_horizon(ii,jj), jj=1,ncols )
     end do
     close(50)
980  format( 10000000(i6, 1x) )
     !
     print*, 'processing soil-id for horizon completed', hh
  end do
  deallocate(soilid_horizon)
  deallocate(soilid)

  !!
  !! now start consolidating info. for the bulk density stats
  do jj = 1, ncols_out
     do ii = 1, nrows_out
        if( grid_out(ii,jj)%nsoils < 1 ) cycle
        !! allocate space for consolidating individual bulk dens.
        allocate( grid_out(ii,jj)%bd(grid_out(ii,jj)%nsoils) )
        grid_out(ii,jj)%bd(:) = 1.5_sp
     end do
  end do
  print*, 'allocating size for consolidation over'
  !!
  !! once again go through soil LUT consolidating bulk density data
  grid_out(:,:)%nsoils = 0  !! use this variable again as a counter
  do kk = 1, nsoil_types
      do hh = 1, nsoil_layers
        sand         = soil_orig(kk)%sand(hh)
        clay         = soil_orig(kk)%clay(hh) 
        bulk_density = soil_orig(kk)%bd  (hh) 
        grid_out(sand,clay)%nsoils = grid_out(sand,clay)%nsoils + 1
        grid_out(sand,clay)%bd(grid_out(sand,clay)%nsoils) = bulk_density
     end do
  end do
  print*, 'consolidating LUT over'
  !!
  !!  calculate mean and standard deviation
  grid_out(:,:)%bd_mean  = -9999.0_sp
  grid_out(:,:)%bd_stdev = -9999.0_sp
  do jj = 1, ncols_out
     do ii = 1, nrows_out
        if( grid_out(ii,jj)%nsoils == 0) cycle
        if( grid_out(ii,jj)%nsoils < 2 ) then
           grid_out(ii,jj)%bd_mean  = grid_out(ii,jj)%bd(1) 
           grid_out(ii,jj)%bd_stdev = 0.001_sp
        else
           grid_out(ii,jj)%bd_mean  = average( grid_out(ii,jj)%bd(:) )
           grid_out(ii,jj)%bd_stdev =  stddev( grid_out(ii,jj)%bd(:) )
        end if
        !! deallocate space
        deallocate( grid_out(ii,jj)%bd )
        !!
     end do
  end do
  print*, 'mean and std dev. of Bd over'
  !!
  !! clean up data outside of mask
  where( grid_out%soil_id < 1)
     grid_out%soil_id  = int(nodata_value_out)
     grid_out%nsoils   = int(nodata_value_out)
     grid_out%bd_mean  = nodata_value_out
     grid_out%bd_stdev = nodata_value_out
  end where
  !!
  !! WRITE ASCII FILES
  fName = trim(dataPath)//trim('soil_id_reduced_grid.asc') 
  open (50, file= trim(fName), status='unknown', action = 'write')
  write(50,1)'ncols       ', ncols_out
  write(50,1)'nrows       ', nrows_out
  write(50,2)'xllcorner   ', xllcorner_out
  write(50,2)'yllcorner   ', yllcorner_out
  write(50,2)'cellsize    ', cellsize_out
  write(50,1)'NODATA_value', int(nodata_value_out)
  do ii=1,nrows_out
     write(50, 981) ( grid_out(ii,jj)%soil_id, jj=1,ncols_out )
  end do
  close(50)
  !!
  !! write
  fName = trim(dataPath)//trim('soil_frequency_reduced_grid.asc') 
  open (50, file= trim(fName), status='unknown', action = 'write')
  write(50,1)'ncols       ', ncols_out
  write(50,1)'nrows       ', nrows_out
  write(50,2)'xllcorner   ', xllcorner_out
  write(50,2)'yllcorner   ', yllcorner_out
  write(50,2)'cellsize    ', cellsize_out
  write(50,1)'NODATA_value', int(nodata_value_out)
  do ii=1,nrows_out
     write(50, 981) ( grid_out(ii,jj)%nsoils, jj=1,ncols_out )
  end do
  close(50)
  !!
  !! write
  fName = trim(dataPath)//trim('bd_mean_reduced_grid.asc') 
  open (50, file= trim(fName), status='unknown', action = 'write')
  write(50,1)'ncols       ', ncols_out
  write(50,1)'nrows       ', nrows_out
  write(50,2)'xllcorner   ', xllcorner_out
  write(50,2)'yllcorner   ', yllcorner_out
  write(50,2)'cellsize    ', cellsize_out
  write(50,1)'NODATA_value', nodata_value_out
  do ii=1,nrows_out
     write(50, 982) ( grid_out(ii,jj)%bd_mean, jj=1,ncols_out )
  end do
  close(50)
  !!
  !! WRITE  
  fName = trim(dataPath)//trim('bd_stdev_reduced_grid.asc') 
  open (50, file= trim(fName), status='unknown', action = 'write')
  write(50,1)'ncols       ', ncols_out
  write(50,1)'nrows       ', nrows_out
  write(50,2)'xllcorner   ', xllcorner_out
  write(50,2)'yllcorner   ', yllcorner_out
  write(50,2)'cellsize    ', cellsize_out
  write(50,1)'NODATA_value', nodata_value_out
  do ii=1,nrows_out
     write(50, 982) ( grid_out(ii,jj)%bd_stdev, jj=1,ncols_out )
  end do
  close(50)
  !!
  !! formats
1 format( a12, 2x, i15   )
2 format( a12, 2x, f15.6 )
981 format( 100(i10,   1x) )
982 format( 100(f10.3, 1x) )
  !!
  !! write results -- LUT
  fName = trim(dataPath)//trim('soil_classdefinition_reduced_grid.txt') 
  open(40, file = trim(fName), status = 'unknown', action = 'write')
  write(40,400)'nSoil_Types', nsoil_types_reduced
  write(40,401)'ID', 'CLAY[%]', 'SAND[%]', 'Bd_mu[gcm-3]', 'Bd_sd[gcm-3]', 'no.Sample'
  !
  do jj = 1, ncols_out
     do ii = 1, nrows_out
        if( grid_out(ii,jj)%soil_id == int(nodata_value_out) ) cycle
        write(40, 402) grid_out(ii,jj)%soil_id, jj, ii, grid_out(ii,jj)%bd_mean, grid_out(ii,jj)%bd_stdev, grid_out(ii,jj)%nsoils
     end do
  end do
  close(40)
  !!
  !! formats
  400 format(A15, I15)
  401 format(A10, 2A12, 2A15  , A10)
  402 format(I10, 2I12, 2F15.3, I10)


  
  !!
 end program main




