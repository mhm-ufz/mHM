! ------------------------------------------------------------------------------
!
! Main Program mtclim preprocessor
!
! author: Johannes Brenner
!
! created: 15.06.2017
!
! ------------------------------------------------------------------------------
!
program main
!
use mo_kind,        only: i4, dp
use mo_mtclim
use mo_time,        only: d2o, jd, fmtdate, get_time
use mo_netcdf
use mo_read_spatial_data
use mo_ncread,      only: Get_NcVarAtt
!
!----------------------------------------
! variable declarations namelists
!
real(dp), dimension(14) :: parameterset
! ini namefile
! input, output file names
character(256) :: in_latlon, in_pre, in_tmax, in_tmin, in_dem, in_asp, in_slp, &
   in_alb_bs, in_alb_ws, out
namelist / io_ini / in_latlon, in_pre, in_tmax, in_tmin, in_dem, in_asp, in_slp, &
   in_alb_bs, in_alb_ws, out
!
! controls
integer(i4) :: outhum, netrad, lwrad, netlwrad
!
namelist / control_ini / outhum, netrad, lwrad, netlwrad
!
! site characteristics
!real(dp)              :: base_elev,base_isoh,site_lat,site_elev,site_slp,site_asp,site_isoh
real(dp)              :: site_ehoriz,site_whoriz,tmax_lr,tmin_lr
!
namelist / parameters_ini / site_ehoriz,site_whoriz,tmax_lr,tmin_lr
!base_elev,base_isoh,site_lat,site_elev,site_slp,site_asp,site_isoh
! from external data
!
! parameters
!
real(dp)      :: TBASE,ABASE,C,B0,B1,B2,RAIN_SCALAR,DIF_ALB,SC_INT,SC_SLOPE, &
SNOW_TCRIT,SNOW_TRATE,TDAYCOEF,LWAVE_COR
namelist / parameters / TBASE,ABASE,C,B0,B1,B2,RAIN_SCALAR,DIF_ALB,SC_INT,SC_SLOPE, &
     SNOW_TCRIT,SNOW_TRATE,TDAYCOEF,LWAVE_COR
!
!----------------------------------------
! variable declarations
!
! netCDF read
!real(dp),     dimension(:,:,:), allocatable :: nc_tmax
!real(dp),     dimension(:,:,:), allocatable :: nc_tmin
!real(dp),     dimension(:,:,:), allocatable :: nc_prcp
! read .asc
integer(i4)                 :: header_nRows, header_nCols, xres, yres
integer(i4)                 :: yres_alb, xres_alb
integer(i4)                 :: month
real(dp)                    :: header_cellsize, header_xllcorner, header_yllcorner
real(dp)                    :: header_nodata
real(dp),    dimension(:,:), allocatable :: data_dem, data_slp, data_asp
logical,     dimension(:,:), allocatable :: mask_asc
real(dp),    dimension(:,:), allocatable :: dem_cell, asp_cell, slp_cell
real(dp),    dimension(:,:), allocatable :: alb_ws_cell, alb_bs_cell
real(dp)                                 :: dem, asp, slp

! netCDF write
type(NcDataset)          :: nc
type(NcDimension)        :: dim_y, dim_x, dim_time
type(NcVariable)         :: var_lat, var_lon, var_x, var_y, var_time
type(NcVariable)         :: var_pre, var_tmin, var_tmax
type(NcVariable)         :: var_eabs, var_swin, var_swrad_dif
type(NcVariable)         :: var_net_rad, var_sradnet, var_lwrad, var_lradnet
type(NcVariable)         :: var_albedo
type(NcVariable)         :: var_alb_ws, var_alb_bs
!
character(*), parameter  :: vname_time="time", vname_lat="lat", vname_lon="lon"
character(*), parameter  :: vname_y="yc", vname_x="xc"
character(*), parameter  :: vname_pre="pre", vname_tmin="tmin", vname_tmax="tmax"
character(*), parameter  :: vname_eabs="eabs", vname_vpd="vpd"
character(*), parameter  :: vname_swin="swrad_in", vname_swdif="swrad_diffuse"
character(*), parameter  :: vname_net_rad="net_rad", vname_sradnet="net_swrad"
character(*), parameter  :: vname_lwrad="lwrad_in", vname_lradnet="net_lwrad"
character(*), parameter  :: vname_alb_ws="albedo_ws", vname_alb_bs="albedo_bs", vname_albedo="albedo"
real(dp),    allocatable :: wlat(:,:), wlon(:,:)
integer(i4), allocatable :: wtime(:), wxc(:), wyc(:)
real(dp)                 :: wfvalue
real(dp),    allocatable :: nc_tmax(:,:,:), nc_tmin(:,:,:), nc_prcp(:,:,:)
real(dp),    allocatable :: nc_alb_bs(:,:,:), nc_alb_ws(:,:,:)

! datetime
integer(i4)                             :: julStart, julEnd
real(dp)                                :: jul_dp
integer(i4)                             :: date(8)
character(i4)                           :: timestr(1)
character(256)                          :: timeunit
logical,      dimension(:), allocatable :: mask
real(dp)                                :: albedo
! mtCLIM
integer(i4)                             :: ndays, ndays_m
integer(i4),  dimension(:), allocatable :: yday         !/* array of yearday values */
integer(i4),  dimension(:), allocatable :: yday_m       !/* array of yearday values / masked */
real(dp),     dimension(:), allocatable :: tmax_f       !/* array of base maximum temperature values */
real(dp),     dimension(:), allocatable :: tmin_f       !/* array of base minimum temperature values */
real(dp),     dimension(:), allocatable :: prcp_f       !/* array of base daily precipitation values */
real(dp),     dimension(:), allocatable :: tmax         !/* array of base maximum temperature values / masked */
real(dp),     dimension(:), allocatable :: tmin         !/* array of base minimum temperature values / masked */
real(dp),     dimension(:), allocatable :: prcp         !/* array of base daily precipitation values / masked */
real(dp),     dimension(:), allocatable :: s_tmax       !/* array of site tmax values */
real(dp),     dimension(:), allocatable :: s_tmin       !/* array of site tmin values */
real(dp),     dimension(:), allocatable :: s_tday       !/* array of site daylight temperature values */
real(dp),     dimension(:), allocatable :: s_prcp       !/* array of site prcp values */
real(dp),     dimension(:), allocatable :: s_hum        !/* array of site humidity values (VPD or VP, Pa) */
real(dp),     dimension(:), allocatable :: s_srad       !/* array of site shortwave radiation values */
real(dp),     dimension(:), allocatable :: s_srad_dif   !/* array of site shortwave radiation values */
real(dp),     dimension(:), allocatable :: s_lrad       !/* array of site longwave radiation values */
real(dp),     dimension(:), allocatable :: s_sradnet    !/* array of site net shortwave radiation values */
real(dp),     dimension(:), allocatable :: s_lradnet    !/* array of site net longwave radiation values */
real(dp),     dimension(:), allocatable :: s_dayl       !/* array of site daylength values */
real(dp),     dimension(:), allocatable :: s_swe        !/* array of site snow water equivalent values (cm) */
!
! output
real(dp),     dimension(:,:,:), allocatable :: out_eabs     !/* array of humidity values (VP - Pa) */
real(dp),     dimension(:,:,:), allocatable :: out_swrad    !/* array of shortwave radiation (W m-2) */
real(dp),     dimension(:,:,:), allocatable :: out_srad_dif !/* array of shortwave radiation (W m-2) */
!
! additional optional output
real(dp),     dimension(:,:,:), allocatable :: out_net_rad  !/* array of site net radiation values (W m-2) */
real(dp),     dimension(:,:,:), allocatable :: out_sradnet  !/* array of site net longwave radiation values (W m-2) */
real(dp),     dimension(:,:,:), allocatable :: out_lradnet  !/* array of site net longwave radiation values (W m-2) */
real(dp),     dimension(:,:,:), allocatable :: out_lwrad    !/* array of longwave radiation (W m-2) */
real(dp),     dimension(:,:,:), allocatable :: out_albedo   !/* array of albedo (--) */
!
integer(i4),  dimension(3)              :: dim_size
!
integer(i4)                             :: ii,row,col
!
logical     :: isgood, allgood
allgood = .true.
!
Write(*,*) ''
Write(*,*) 'preprocessor mo_mtclim.f90'
!
!----------------------------------------
! read out namefiles ini, parameters
!
open(700, file='ini')
read(700, nml=io_ini)
read(700, nml=control_ini)
read(700, nml=parameters_ini)
close(700)
!
open(700, file='parameters')
read(700, nml=parameters)
close(700)
!
parameterset(1) = TBASE
parameterset(2) = ABASE
parameterset(3) = C
parameterset(4) = B0
parameterset(5) = B1
parameterset(6) = B2
parameterset(7) = RAIN_SCALAR
parameterset(8) = DIF_ALB
parameterset(9) = SC_INT
parameterset(10) = SC_SLOPE
parameterset(11) = SNOW_TCRIT
parameterset(12) = SNOW_TRATE
parameterset(13) = TDAYCOEF
parameterset(14) = LWAVE_COR
!
!----------------------------------------
! read nc files pre.nc, tmin.nc, tmax.nc
!
!----------------------------------------
! precipitation
! read time, lat/lon and data
! open dataset
nc = NcDataset(in_latlon,"r")
! access the variable
var_lat  = nc%getVariable(vname_lat)
var_lon  = nc%getVariable(vname_lon)
var_y = nc%getVariable(vname_y)
var_x = nc%getVariable(vname_x)
! read the data
call var_lat%getData(wlat)
call var_lon%getData(wlon)
call var_y%getData(wyc)
call var_x%getData(wxc)
! close dataset
call nc%close()
! pre
nc = NcDataset(in_pre,"r")
! access the variable
var_time = nc%getVariable(vname_time)
var_pre  = nc%getVariable(vname_pre)
! read the data
call var_time%getData(wtime)
call var_pre%getData(nc_prcp)
! read the fill value
call var_pre%getFillValue(wfvalue)
! close dataset
call nc%close()
!
! get time unit argument
call Get_NcVarAtt(in_pre, 'time', 'units', timeunit)

!----------------------------------------
! maximum air temperatures
! open dataset
nc = NcDataset(in_tmax,"r")
! access the variable
var_tmax = nc%getVariable(vname_tmax)
! read the data
call var_tmax%getData(nc_tmax)
! close dataset
call nc%close()
!----------------------------------------
! minimum air temperatures
! open dataset
nc = NcDataset(in_tmin,"r")
! access the variable
var_tmin = nc%getVariable(vname_tmin)
! read the data
call var_tmin%getData(nc_tmin)
! close dataset
call nc%close()
!----------------------------------------
! get dimension
!dim = Get_NcDim(Filename=in_pre, Variable='pre')
! x dimension - lon - col - x
!dim_size(1) = size(wxc)
dim_size(1) = size(nc_prcp,1)
! y dimension - lat - row - y
!dim_size(2) = size(wyc)
dim_size(2) = size(nc_prcp,2)
! time dimension
dim_size(3) = size(wtime)
ndays = dim_size(3)
! get datetime julian
call get_time(in_pre, 'pre', julStart, julEnd)
! create julian time series
allocate(yday(ndays))
do ii = 0, ndays-1, 1
  jul_dp = julStart + ii
  date = jd(jul_dp)
  timestr = fmtdate(date,'%O')
  read(timestr(1),'(I3)') yday(ii+1)
end do
!----------------------------------------
! read .asc files dem.asc, aspect.asc, slope.asc
!
call read_header_ascii(in_dem,777, &
   header_nCols, header_nRows,  &
   header_xllcorner, header_yllcorner, &
   header_cellsize, header_nodata)
! dem read
call read_spatial_data_ascii(in_dem, 778, &
   header_nCols, header_nRows,  &
   header_xllcorner, header_yllcorner, header_cellsize, &
   data_dem, mask_asc)
! aspect read
call read_spatial_data_ascii(in_asp, 779, &
   header_nCols, header_nRows, &
   header_xllcorner, header_yllcorner, header_cellsize, &
   data_asp, mask_asc)
! slope read
call read_spatial_data_ascii(in_slp, 780, &
   header_nCols, header_nRows, &
   header_xllcorner, header_yllcorner, header_cellsize, &
   data_slp, mask_asc)

!----------------------------------------
! albedo netCDF (white-sky - ws, black-sky - bs)
if ( netrad .gt. 0_i4 ) then
  ! albedo from 12moth climatology
  ! read in albedo netCDF
  !
  ! white-sky albedo
  ! open dataset
  nc = NcDataset(in_alb_ws,"r")
  ! access the variable
  var_alb_ws = nc%getVariable(vname_alb_ws)
  ! read the data
  call var_alb_ws%getData(nc_alb_ws)
  ! close dataset
  call nc%close()
  !
  ! black-sky albedo
  ! open dataset
  nc = NcDataset(in_alb_bs,"r")
  ! access the variable
  var_alb_bs = nc%getVariable(vname_alb_bs)
  ! read the data
  call var_alb_bs%getData(nc_alb_bs)
  ! close dataset
  call nc%close()
  !
  ! get res x,y of albeo data
  yres_alb = size(nc_alb_bs, dim=2) / dim_size(2)
  xres_alb = size(nc_alb_bs, dim=1) / dim_size(1)
  !
  ! allocate space
  allocate(alb_ws_cell(yres_alb,xres_alb))
  allocate(alb_bs_cell(yres_alb,xres_alb))
  !
end if
!----------------------------------------
! old readin nc
! allocate nc variables
!allocate(nc_prcp(dim(1),dim(2),dim(3)))
!allocate(nc_tmin(dim(1),dim(2),dim(3)))
!allocate(nc_tmax(dim(1),dim(2),dim(3)))
! read in data
!call Get_NcVar(Filename=in_pre, VarName='pre', Dat=nc_prcp)
!call Get_NcVar(Filename=in_tmin, VarName='tmin', Dat=nc_tmin)
!call Get_NcVar(Filename=in_tmax, VarName='tmax', Dat=nc_tmax)
!
!print *, nc_prcp(100,100,100)
!print *, nc_tmin(100,100,100)
!print *, nc_tmax(100,100,100)
!
!----------------------------------------
! allocate space in the data arrays for input data
allocate(mask(ndays))
allocate(tmax_f(ndays))
allocate(tmin_f(ndays))
allocate(prcp_f(ndays))
!----------------------------------------
! allocate space in out arrays
allocate(out_eabs(dim_size(1),dim_size(2),ndays))
allocate(out_swrad(dim_size(1),dim_size(2),ndays))
allocate(out_srad_dif(dim_size(1),dim_size(2),ndays))
if (netrad .gt. 0_i4) then
  allocate(out_net_rad(dim_size(1),dim_size(2),ndays))
  allocate(out_sradnet(dim_size(1),dim_size(2),ndays))
  allocate(out_albedo(dim_size(1),dim_size(2),ndays))
end if
if (netlwrad .gt. 0_i4) then
  allocate(out_lradnet(dim_size(1),dim_size(2),ndays))
end if
if (lwrad .gt. 0_i4) then
  allocate(out_lwrad(dim_size(1),dim_size(2),ndays))
end if
!
!----------------------------------------
! loop over meteo field time series
!
!
yres = header_nRows / dim_size(2)
xres = header_nCols / dim_size(1)
!
allocate(dem_cell(xres,yres))
allocate(asp_cell(xres,yres))
allocate(slp_cell(xres,yres))
!
do row = 1, dim_size(1), 1
  do col = 1, dim_size(2), 1
    !
    ! aggregate .asc data
    dem_cell = data_dem((1+(row-1)*yres):row*yres, (1+(col-1)*xres):col*xres)
    asp_cell = data_asp((1+(row-1)*yres):row*yres, (1+(col-1)*xres):col*xres)
    slp_cell = data_slp((1+(row-1)*yres):row*yres, (1+(col-1)*xres):col*xres)
    if (all(dem_cell .le. wfvalue)) then
    !if (any(dem_cell .le. wfvalue)) then
      dem = wfvalue
      slp = wfvalue
      asp = wfvalue
    else
      dem = sum(dem_cell, dem_cell .ge. 0_dp)/(max(1,size(pack(dem_cell, dem_cell .ge. 0_dp))))
      slp = sum(slp_cell, slp_cell .ge. 0_dp)/(max(1,size(pack(slp_cell, dem_cell .ge. 0_dp))))
      asp = sum(asp_cell, asp_cell .ge. 0_dp)/(max(1,size(pack(asp_cell, dem_cell .ge. 0_dp))))
    end if
    ! fill full meteo 1D array
    prcp_f = nc_prcp(row,col,:)
    tmin_f = nc_tmin(row,col,:)
    tmax_f = nc_tmax(row,col,:)
    ! NaN masking, all variables
    do ii = 1, ndays, 1
      mask(ii) = (prcp_f(ii) .le. wfvalue) .or. &
      (tmin_f(ii) .le. wfvalue) .or. (tmax_f(ii) .le. wfvalue)
    end do
    ! if no valid data point in 1D array
    ! write NaN - -9999._dp
    if (all(mask) .or. (dem .le. wfvalue)) then
      !
      print *, 'no valid data point - nothing to do'
      !
    else
      ! reduce length of input
      ! update ndays
      ndays_m = size(pack(prcp_f, .not. mask))
      ! allocate masked size
      allocate(tmin(ndays_m))
      allocate(tmax(ndays_m))
      allocate(prcp(ndays_m))
      allocate(yday_m(ndays_m))
      !
      allocate(s_swe(ndays_m))
      !
      allocate(s_tmin(ndays_m))
      allocate(s_tmax(ndays_m))
      allocate(s_tday(ndays_m))
      !
      allocate(s_hum(ndays_m))
      allocate(s_srad(ndays_m))
      allocate(s_srad_dif(ndays_m))
      allocate(s_lrad(ndays_m))
      allocate(s_lradnet(ndays_m))
      allocate(s_sradnet(ndays_m))
      allocate(s_dayl(ndays_m))
      !
      ! fill variables
      do ii = 1, ndays_m, 1
        tmin = pack(tmin_f,.not. mask)
        tmax = pack(tmax_f,.not. mask)
        prcp = pack(prcp_f,.not. mask)
        yday_m = pack(yday,.not. mask)
      end do
      !
      ! run mtCLIM
      call eval_mtclim (parameterset, yday_m, tmin, tmax, prcp, ndays_m, &
           outhum, lwrad, netlwrad, &
           s_tmax, s_tmin, s_tday, s_prcp, &
           s_hum, s_srad, s_srad_dif, s_lrad, s_lradnet, s_dayl, &
           site_elev=dem, base_elev=dem, tmin_lr=tmin_lr, tmax_lr=tmax_lr, &
           site_isoh=50._dp, base_isoh=50._dp, &
           site_lat=wlat(row,col), site_asp=asp, site_slp=slp, &
           site_ehoriz=site_ehoriz, site_whoriz=site_whoriz)
      !
      deallocate(tmin)
      deallocate(tmax)
      deallocate(prcp)
      deallocate(yday_m)
      deallocate(s_swe)
    end if
    !
    ! net shortwave radiation - albedo
    if (all(mask) .or. (dem .le. wfvalue)) then
      if (netrad .gt. 0_i4) then
        out_net_rad(row,col,:) = wfvalue
        out_sradnet(row,col,:) = wfvalue
        out_albedo(row,col,:)  = wfvalue
      end if
      if (netlwrad .gt. 0_i4) then
        out_lradnet(row,col,:) = wfvalue
      end if
      out_eabs(row,col,:) = wfvalue
      out_swrad(row,col,:) = wfvalue
      out_srad_dif(row,col,:) = wfvalue
      if (lwrad .gt. 0_i4) then
        out_lwrad(row,col,:) = wfvalue
      end if
    else
      ! write to out object
      out_eabs(row,col,:) = s_hum(:)
      out_swrad(row,col,:) = s_srad(:)
      out_srad_dif(row,col,:) = s_srad_dif(:)
      ! if albedo is given calc net radiation
      if (netrad .gt. 0_i4) then
        ! upscale albedo data to meteo resolution
        !
        ! calc month of current date
        do ii = 1, ndays_m, 1
          jul_dp = julStart + ii - 1_dp
          date = jd(jul_dp)
          month = date(2)
          ! select case month
          ! and aggregate albebo data
          !
          select case (month)
          case(1)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 1)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 1)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(2)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 2)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 2)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(3)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 3)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 3)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(4)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 4)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 4)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(5)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 5)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 5)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(6)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 6)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 6)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(7)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 7)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 7)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(8)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 8)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 8)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(9)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 9)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 9)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(10)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 10)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 10)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(11)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 11)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 11)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          case(12)
            alb_ws_cell = nc_alb_ws((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 12)
            alb_bs_cell = nc_alb_bs((1+(row-1)*yres_alb):row*yres_alb, (1+(col-1)*xres):col*xres_alb, 12)
            if (all(alb_ws_cell .le. wfvalue) .or. all(alb_bs_cell .le. wfvalue)) then
            !if (any(alb_ws_cell .le. wfvalue) .or. any(alb_bs_cell .le. wfvalue)) then
              alb_ws = wfvalue
              alb_bs = wfvalue
            else
              alb_ws = sum(alb_ws_cell, alb_ws_cell .ge. 0_dp)/(max(1,size(pack(alb_ws_cell, dem_cell .ge. 0_dp))))
              alb_bs = sum(alb_bs_cell, alb_bs_cell .ge. 0_dp)/(max(1,size(pack(alb_bs_cell, dem_cell .ge. 0_dp))))
            end if
          end select
          !
          ! calc albedo in current month for net rad calc
          ! [Lewis and Barnsley,1994]: α(Ωi) = α_ws(Ωi) * p + α_bs(Ωi) * (1-p)
          ! p: ratio of the surface downward diffuse shortwave radiation to
          !    the surface downward total shortwave radiation
          ! This simple formula assumes the isotropic reflectance of the surface
          ! and fails to incorporate multiple scattering between the surface
          ! and the atmosphere [Román et al., 2010].
          ! Román et al. [2010] developed a group of new equations to address these issues.
          !
          albedo = alb_ws * s_srad_dif(ii) / s_srad(ii) + alb_bs * &
          ( 1.0_dp - s_srad_dif(ii) / s_srad(ii))
          out_albedo(row,col,ii)  = albedo
          out_sradnet(row,col,ii) = s_srad(ii) *albedo
          out_net_rad(row,col,ii) = s_srad(ii) *albedo + s_lradnet(ii)
        end do
      end if
      ! if net longwave radiation is calculated return it
      if (netlwrad .gt. 0_i4) then
        out_lradnet(row,col,:) = s_lradnet(:)
      end if
      ! if longwave radiation is calculated return it
      if (lwrad .gt. 0_i4) then
        out_lwrad(row,col,:) = s_lrad(:)
      end if
      !
    end if
    !
    if (.not. all(mask) .and. (dem .gt. wfvalue)) then
      !deallocate mtclim variables
      deallocate(s_tmin)
      deallocate(s_tmax)
      deallocate(s_tday)
      deallocate(s_hum)
      deallocate(s_srad)
      deallocate(s_srad_dif)
      deallocate(s_lrad)
      deallocate(s_lradnet)
      deallocate(s_sradnet)
      deallocate(s_dayl)
      !
    end if
  end do
end do
!----------------------------------------
! write output to netCDF file
!
! create a file
nc = NcDataset(out, "w")
! create dimensions
dim_x    = nc%setDimension("x", dim_size(1))
dim_y    = nc%setDimension("y", dim_size(2))
dim_time = nc%setDimension("time", -1) ! length < 0 -> unlimited dimension
! create variables
var_time = nc%setVariable(vname_time, "i32", (/dim_time/))
!var_lat  = nc%setVariable(vname_lat,  "f32", (/dim_x, dim_y/))
!var_lon  = nc%setVariable(vname_lon , "f32", (/dim_x, dim_y/))
var_x  = nc%setVariable("x" , "f32", (/dim_x/))
var_y  = nc%setVariable("y",  "f32", (/dim_y/))

! actual vapore pressure (Pa) or vapore pressure deficit (Pa)
select case (outhum)
case(0)
  var_eabs = nc%setVariable(vname_vpd, "f64", (/dim_x, dim_y, dim_time/))
case(1)
  var_eabs = nc%setVariable(vname_eabs, "f64", (/dim_x, dim_y, dim_time/))
end select
! incoming shortwave radiation (W/m2)
var_swin = nc%setVariable(vname_swin, "f64", (/dim_x, dim_y, dim_time/))
! diffuse shortwave radiation (W/m2)
var_swrad_dif = nc%setVariable(vname_swdif, "f64", (/dim_x, dim_y, dim_time/))
!
if (lwrad .gt. 0_i4) then
  ! longwave radiation
  var_lwrad = nc%setVariable(vname_lwrad, "f64", (/dim_x, dim_y, dim_time/))
end if
!
if (netrad .gt. 0_i4) then
  ! net radiation (W/m2)
  var_net_rad = nc%setVariable(vname_net_rad, "f64", (/dim_x, dim_y, dim_time/))
  ! albedo (--)
  var_albedo = nc%setVariable(vname_albedo, "f64", (/dim_x, dim_y, dim_time/))
  ! net shortwave radiation
  var_sradnet = nc%setVariable(vname_sradnet, "f64", (/dim_x, dim_y, dim_time/))
end if
!
if (netlwrad .gt. 0_i4) then
  ! net longwave radiation
  var_lradnet = nc%setVariable(vname_lradnet, "f64", (/dim_x, dim_y, dim_time/))
end if
!
! add some variable attributes
call var_time%setAttribute("units", timeunit)
! set fill value before any data is written
call var_eabs%setFillValue(wfvalue)
call var_swin%setFillValue(wfvalue)
call var_swrad_dif%setFillValue(wfvalue)
if (netrad .gt. 0_i4) then
  call var_net_rad%setFillValue(wfvalue)
  call var_albedo%setFillValue(wfvalue)
  call var_sradnet%setFillValue(wfvalue)
end if
if (netlwrad .gt. 0_i4) then
  call var_lradnet%setFillValue(wfvalue)
end if
if (lwrad .gt. 0_i4) then
  call var_lwrad%setFillValue(wfvalue)
end if
! write data of static variables
!call var_lat%setData(wlat)
!call var_lon%setData(wlon)
call var_y%setData(wyc)
call var_x%setData(wxc)
! append data within a loop
do ii=1, ndays_m, 1
  call var_time%setData(wtime(ii), start=(/ii/))
  ! actual vapore pressure (Pa)
  call var_eabs%setData(out_eabs(:,:,ii), start=(/1,1,ii/))
  ! incoming shortwave radiation (W/m2)
  call var_swin%setData(out_swrad(:,:,ii), start=(/1,1,ii/))
  ! diffuse shortwave radiation (W/m2)
  call var_swrad_dif%setData(out_srad_dif(:,:,ii), start=(/1,1,ii/))
  !
  if (netrad .gt. 0_i4) then
    ! net radiation (W/m2)
    call var_net_rad%setData(out_net_rad(:,:,ii), start=(/1,1,ii/))
    call var_albedo%setData(out_albedo(:,:,ii),   start=(/1,1,ii/))
    call var_sradnet%setData(out_sradnet(:,:,ii), start=(/1,1,ii/))
  end if
  !
  if (netlwrad .gt. 0_i4) then
    ! net longwave radiation (W/m2)
    call var_lradnet%setData(out_lradnet(:,:,ii), start=(/1,1,ii/))
  end if
  ! longwave radiation (w/m2)
  if (lwrad .gt. 0_i4) then
    call var_lwrad%setData(out_lwrad(:,:,ii), start=(/1,1,ii/))
  end if
end do

! add some more variable attributes
! data
! actual vapore pressure (Pa)
select case (outhum)
case(0)
  call var_eabs%setAttribute("standard_name", vname_vpd)
  call var_eabs%setAttribute("long_name",   "Vapore pressure deficit")
case(1)
  call var_eabs%setAttribute("standard_name", vname_eabs)
  call var_eabs%setAttribute("long_name",   "Vapore pressure actual")
end select
call var_eabs%setAttribute("units",   "Pa")
call var_eabs%setAttribute("missing_value", wfvalue)
! incoming shortwave radiation (W/m2)
call var_swin%setAttribute("standard_name", vname_swin)
call var_swin%setAttribute("long_name",   "Incoming shortwave radiation")
call var_swin%setAttribute("units",   "W m-1")
call var_swin%setAttribute("missing_value", wfvalue)
! diffuse shortwave radiation (W/m2)
call var_swrad_dif%setAttribute("standard_name", vname_swdif)
call var_swrad_dif%setAttribute("long_name",   "Diffuse shortwave radiation")
call var_swrad_dif%setAttribute("units",   "W m-1")
call var_swrad_dif%setAttribute("missing_value", wfvalue)
!
if (netrad .gt. 0_i4) then
  ! net radiation (W/m2)
  call var_net_rad%setAttribute("standard_name", vname_net_rad)
  call var_net_rad%setAttribute("long_name",   "Net radiation")
  call var_net_rad%setAttribute("units",   "W m-2")
  call var_net_rad%setAttribute("missing_value", wfvalue)
  ! albedo (--)
  call var_albedo%setAttribute("standard_name", vname_albedo)
  call var_albedo%setAttribute("long_name",   "Albedo")
  call var_albedo%setAttribute("units",   "-")
  call var_albedo%setAttribute("missing_value", wfvalue)
  ! net shortwave radiation (--)
  call var_sradnet%setAttribute("standard_name", vname_sradnet)
  call var_sradnet%setAttribute("long_name",   "Net shortwave radiation")
  call var_sradnet%setAttribute("units",   "W m-2")
  call var_sradnet%setAttribute("missing_value", wfvalue)
end if
!
if (netlwrad .gt. 0_i4) then
  ! net longwave radiation (W/m2)
  call var_lradnet%setAttribute("standard_name", vname_lradnet)
  call var_lradnet%setAttribute("long_name",   "Net longwave radiation")
  call var_lradnet%setAttribute("units",   "W m-2")
  call var_lradnet%setAttribute("missing_value", wfvalue)
end if
!
if (lwrad .gt. 0_i4) then
  ! longwave radiation (W/m2)
  call var_lwrad%setAttribute("standard_name", vname_lwrad)
  call var_lwrad%setAttribute("long_name",   "Longwave radiation")
  call var_lwrad%setAttribute("units",   "W m-2")
  call var_lwrad%setAttribute("missing_value", wfvalue)
end if
! time
call var_time%setAttribute("standard_name", "time")
call var_time%setAttribute("long_name", "Time variable")
call var_time%setAttribute("calendar", "standard")
call var_time%setAttribute("axis", "T")
! x
call var_x%setAttribute("standard_name", "x coordinate of projection")
call var_x%setAttribute("long_name", "x coordinate of projection")
call var_x%setAttribute("units", "Meter")
call var_x%setAttribute("axis", "X")
! y
call var_y%setAttribute("standard_name", "y coordinate of projection")
call var_y%setAttribute("long_name", "y coordinate of projection")
call var_y%setAttribute("units", "Meter")
call var_y%setAttribute("axis", "Y")
! add global attributes
call nc%setAttribute("Author", "Johannes Brenner")
call nc%setAttribute("Institution", "CHS, UFZ Leipzig")
call nc%setAttribute("Title", "retrived from mtCLIM v4.3")
!call nc%setAttribute("creation_date", )

! close the file
call nc%close()

!----------------------------------------
! free previously allocated memory before returning
!
deallocate(yday)
deallocate(tmax_f)
deallocate(tmin_f)
deallocate(prcp_f)
!
if (netrad .gt. 0_i4) then
  deallocate(out_net_rad)
  deallocate(out_sradnet)
end if
!
if (netlwrad .gt. 0_i4) then
  deallocate(out_lradnet)
end if
deallocate(out_eabs)
deallocate(out_swrad)
deallocate(out_srad_dif)
if (lwrad .gt. 0_i4) then
  deallocate(out_lwrad)
end if
!
!------------------------------------------
isgood = .true.
! Finish
allgood = allgood .and. isgood
  if (allgood) then
     write(*,*) 'mo_mtclim o.k.'
  else
     write(*,*) 'mo_mtclim failed!'
  endif
!
end program main
