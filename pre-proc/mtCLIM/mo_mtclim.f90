module mo_mtclim

  use mo_kind,          only: i4, dp
  use mo_constants,     only: LR_STD, T_STD, gravity_dp, R, MA, deg2rad_dp, SECPERRAD, &
       MINDECL, DAYSOFF, RADPERDAY, pi_dp, sigma_dp, P0_dp, cp0_dp, EPS

  implicit none
  !
  public::pulled_boxcar                 ! calculates a moving average of antecedent values in an array
  public::atm_pres                      ! calculates the atmospheric pressure as a function of elevation
  public::snowpack                      ! estimates the accumulation and melt of snow
  public::calc_pet                      ! calculates the potential evapotranspiration for aridity
                                        ! corrections in calc_vpd(), according to Kimball et al., 1997
  public::calc_tair                     ! calculates daily air temperatures
  public::calc_prcp                     ! calculates daily total precipitation
  public::calc_srad_humidity_iterative  ! estimate srad and humidity with iterative algorithm
  public::eval_mtclim                   ! mtclim model evaluation

  !
  !
  ! ------------------------------------------------------------------------------
  !    NAME
  !        mtclim
  !
  !    PURPOSE
  !        uses observations of daily maximum temperature, minimum temperature, and precipitation from
  !        one location (the "base") to estimate the temperature, precipitation, incoming short-wave radiation, and
  !        humidity (vapor pressure deficit) at another location (the "site") or at same location - mtCLIM approach.
  !        The base and the site can be at different elevations, and can have different slopes and aspects.
  !        Additionally, incoming long-wave radiation is estimated using the Prata (1996) clear-sky algorithm,
  !        combined with the Deardorff (1978) full-sky approach (also attributed to Crawford and Duchon, 1999).
  !
  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        see test program in directory test_mo_mtclim

  !    LITERATURE
  !        Thornton, P.E., and S.W. Running, 1999. An improved algorithm for estimating incident daily solar radiation
  !                     from measurements of temperature, humidity, and precipitation.
  !                     Agricultural and Forest Meteorology, 93:211-228.
  !        Bohn, T. J., Livneh, B., Oyler, J. W., Running, S. W., Nijssen, B., & Lettenmaier, D. P. (2013).
  !                     Global evaluation of MTCLIM and related algorithms for forcing of ecological and hydrological models.
  !                     Agricultural and Forest Meteorology, 176, 38–49. http://doi.org/10.1016/j.agrformet.2013.03.003
  !        Deardorff, J. W. (1978). Efficient prediction of ground surface temperature and moisture, with inclusion of
  !                     a layer of vegetation. Journal of Geophysical Research, 83(C4), 1889.
  !                     http://doi.org/10.1029/JC083iC04p01889
  !        Prata, A. J. (1996). A new long-wave formula for estimating downward clear-sky radiation at the surface.
  !                     Quarterly Journal of the Royal Meteorological Society, 122(533), 1127–1151.
  !                     http://doi.org/10.1002/qj.49712253306
  !        Crawford, T. M., Duchon, C. E., Crawford, T. M., & Duchon, C. E. (1999). An Improved Parameterization for Estimating
  !                     Effective Atmospheric Emissivity for Use in Calculating Daytime Downwelling Longwave Radiation.
  !                     Journal of Applied Meteorology, 38(4), 474–480.
  !                     http://doi.org/10.1175/1520-0450(1999)038<0474:AIPFEE>2.0.CO;2

  !    HISTORY
  !        Written, most of the code was written by Peter E. Thornton at the Univeristy of Montana (see below).
  !      	          Adapted by: Johannes Brenner (C to Fortran translation)
  !	                  Adaptation started on Sept 2016 using the mtclim v4.3 code from http://www.ntsg.umt.edu/project/mtclim
  !        Changes, Johannes Brenner, May 2017: longwave radiation calculation according to Bohn et al. 2013, marked as !*/begin longwave, !*/end longwave
  !        Changes, Johannes Brenner, May 2017: constants coming from mo_constants
  !
  !    CHANGES to mtCLIM v4.3
  !        - code includes calculation of incoming longwave radiation (Bohn et al. 2013) and outgoing net longwave radiation (Allen et al. 1999)
  !
  !
  ! ------------------------------------------------------------------------------

  private

  ! ------------------------------------------------------------------------------

contains

  ! ------------------------------------------------------------------------------
  !    NAME
  !        calc_srad_humidity_iterative
  !
  !    PURPOSE
  !        estimate incoming shortwave and longwave radiation,
  !        as well as outgoing net longwave radiation
  !        and humidity/vapore pressure deficit with iterative algorithm,
  !        using maximum temperature, minimum temperature, and precipitation

  !    CALLING SEQUENCE
  !     call calc_srad_humidity_iterative(yday, s_tmax, s_tmin, s_prcp, &
  !                                       s_hum, s_srad, s_lrad, s_lradnet, s_dayl, ndays, outhum, &
  !                                       site_lat, site_elev, site_asp, site_slp, site_ehoriz, site_whoriz, &
  !                                       SECPERRAD, RADPERDEG, RADPERDAY, MINDECL, DAYSOFF, &
  !                                       MA, R, G_STD, T_STD, CP, LR_STD, EPS, PI, SIGMA,&
  !                                       TBASE,ABASE,C,B0,B1,B2,RAIN_SCALAR,DIF_ALB)

  !    INTENT(OUT)
  !        real(dp),     dimension(:) :: s_hum, vapore pressure deficit in

  !    INTENT(OUT)
  !        real(dp),     dimension(:) :: s_srad, incoming shortwave radiation in W/m2, average over daylight period
  !                                              mean s_srad of day+night: s_srad = s_srad * s_dayl / (24*60*60)

  !    INTENT(OUT)
  !        real(dp),     dimension(:) :: s_srad_dif, incoming shortwave diffuse radiation in W/m2, average over daylight period
  !                                                  mean s_srad of day+night: s_srad_dif = s_srad_dif * s_dayl / (24*60*60)

  !    INTENT(OUT)
  !        real(dp),     dimension(:) :: s_lrad, incoming longwave radiation in W/m2

  !    INTENT(OUT)
  !        real(dp),     dimension(:) :: s_lradnet, outgoing net longwave radiation in W/m2

  !    INTENT(OUT)
  !        real(dp),     dimension(:) :: s_dayl, daylength in s (sunrise to sunset, flat horizons)

  !    INTENT(IN)
  !        real(dp)     :: ndays

  !    INTENT(IN)
  !        real(dp),     dimension(:)     :: s_tmax, maximum air temperature in degC (at site)

  !    INTENT(IN)
  !        real(dp),     dimension(:)     :: s_tmin, minimum air temperature in degC (at site)

  !    INTENT(IN)
  !        real(dp),     dimension(:)     :: s_tday, maximum air temperature in degC (at site)

  !    INTENT(IN)
  !        real(dp),     dimension(:)     :: s_prcp, precipitation in mm/day

  !    INTENT(IN)
  !        real(dp),     dimension(:)     :: s_swe, snow water eqivalent in mm

  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        none

  !    LITERATURE
  !        Bristow, K.L., and G.S. Campbell, 1984. On the relationship between incoming solar radiation
  !                      and daily maximum and minimum temperature.
  !                      Agricultural and Forest Meteorology, 31:159-166.
  !        Bohn, T. J., Livneh, B., Oyler, J. W., Running, S. W., Nijssen, B., & Lettenmaier, D. P. (2013).
  !                     Global evaluation of MTCLIM and related algorithms for forcing of ecological and hydrological models.
  !                     Agricultural and Forest Meteorology, 176, 38–49. http://doi.org/10.1016/j.agrformet.2013.03.003
  !        Deardorff, J. W. (1978). Efficient prediction of ground surface temperature and moisture, with inclusion of
  !                     a layer of vegetation. Journal of Geophysical Research, 83(C4), 1889.
  !                     http://doi.org/10.1029/JC083iC04p01889
  !        Prata, A. J. (1996). A new long-wave formula for estimating downward clear-sky radiation at the surface.
  !                     Quarterly Journal of the Royal Meteorological Society, 122(533), 1127–1151.
  !                     http://doi.org/10.1002/qj.49712253306
  !        Crawford, T. M., Duchon, C. E., Crawford, T. M., & Duchon, C. E. (1999). An Improved Parameterization for Estimating
  !                     Effective Atmospheric Emissivity for Use in Calculating Daytime Downwelling Longwave Radiation.
  !                     Journal of Applied Meteorology, 38(4), 474–480.
  !                     http://doi.org/10.1175/1520-0450(1999)038<0474:AIPFEE>2.0.CO;2
  !         Allen, R. G. & Food and Agriculture Organization of the United Nations. (1998).
  !                     Crop evapotranspiration : guidelines for computing crop water requirements.
  !                     Food and Agriculture Organization of the United Nations.
  !                     Retrieved from http://www.fao.org/docrep/X0490E/X0490E00.htm

  !    HISTORY
  !        Written,  Johannes Brenner, Sept 2016
  !        Change,   Added incoming and net longwave radiation computation, Johannes Brenner, Mai 2017
  !        Change,   Constants coming from mo_constants, Johannes Brenner, Mai 2017

  subroutine calc_srad_humidity_iterative( &
       yday, s_tmax, s_tmin, s_tday, s_prcp, s_swe, & ! intent(in), variables
       s_hum, s_srad, s_srad_dif, s_lrad, s_lradnet, s_dayl, &   ! intent(out), variables
       ndays, outhum, lwrad, netlwrad, &
       site_lat, site_elev, site_asp, site_slp, site_ehoriz, site_whoriz, &
       TBASE, ABASE, C, B0, B1, B2, RAIN_SCALAR, DIF_ALB, LWAVE_COR)

  ! input variables
  integer(i4),  dimension(:),  intent(in)  :: yday
  real(dp),     dimension(:),  intent(in)  :: s_tmax, s_tmin, s_tday, s_prcp, s_swe
  integer(i4),  intent(in)                 :: ndays
  ! select case
  integer(i4),  intent(in)                 :: outhum, lwrad, netlwrad
  ! input constants and parameter
  real(dp),     intent(in)                 :: site_lat, site_elev, site_asp
  real(dp),     intent(in)                 :: site_slp, site_ehoriz, site_whoriz
  real(dp),     intent(in)                 :: TBASE, ABASE, C, B0, B1, B2
  real(dp),     intent(in)                 :: RAIN_SCALAR, DIF_ALB, LWAVE_COR
  ! output variables
  real(dp),     dimension(:), intent(out)  :: s_hum, s_srad, s_srad_dif, s_dayl
  real(dp),     dimension(:), intent(out)  :: s_lrad, s_lradnet
  !
  real(dp),     dimension(:), allocatable :: tdew, tmax, tmin
  real(dp),     dimension(:), allocatable :: dtr, sm_dtr
  real(dp),     dimension(:), allocatable :: parray, window, t_fmax
  real(dp),     dimension(:), allocatable :: save_pet
  real(dp),     dimension(:), allocatable :: daylength
  real(dp),     dimension(:), allocatable :: ttmax0
  real(dp),     dimension(:), allocatable :: flat_potrad
  real(dp),     dimension(:), allocatable :: slope_potrad
  real(dp),     dimension(:), allocatable :: s_tskc
  !
  integer(i4)                             :: ami
  !
  real(dp)                                :: xi,E_a  !source emissivity clear sky (Prata approach)
  real(dp)                                :: pva,pvs,vpd
  real(dp)                                :: sum_prcp,ann_prcp,effann_prcp
  real(dp)                                :: t1, t2
  real(dp)                                :: pratio
  real(dp)                                :: lat,coslat,sinlat,dt,h,dh
  real(dp)                                :: cosslp,sinslp,cosasp,sinasp
  real(dp)                                :: decl,cosdecl,sindecl,cosegeom,sinegeom,coshss,hss
  real(dp)                                :: cza,cbsa,coszeh,coszwh
  real(dp)                                :: sc,dir_beam_topa
  real(dp)                                :: sum_flat_potrad, sum_slope_potrad, sum_trans
  real(dp)                                :: cosh_,sinh_
  real(dp)                                :: avg_horizon, slope_excess
  real(dp)                                :: horizon_scalar, slope_scalar
  real(dp)                                :: sky_prop
  real(dp)                                :: t_tmax,b
  real(dp)                                :: t_final,pdif,pdir,srad1,srad2
  real(dp)                                :: pa
  real(dp)                                :: sum_pet,ann_pet
  real(dp)                                :: tmink,pet,ratio,ratio2,ratio3,tdewk
  real(dp)                                :: trans1,trans2,bsg1,bsg2,bsg3,dir_flat_topa,am
  real(dp)                                :: emissivity, cloudfactor
  !
  integer(i4)                             :: start_yday,end_yday, isloop
  integer(i4)                             :: i, j
  ! optical airmass by degrees, defined later
  real(dp),     dimension(:), allocatable :: optam
!
! local array memory allocation
! allocate space for DTR and smoothed DTR arrays
allocate(dtr(ndays))
allocate(sm_dtr(ndays))
! allocate space for effective annual precip array
allocate(parray(ndays))
! allocate space for the prcp totaling array
allocate(window(ndays+90))
! allocate space for t_fmax
allocate(t_fmax(ndays))
! allocate space for save_pet array
allocate(save_pet(ndays))
!allocte space for daylength
allocate(daylength(366))
!allocte space for ttmax0
allocate(ttmax0(366))
!allocte space for flat_potrad
allocate(flat_potrad(366))
!allocte space for slope_potrad
allocate(slope_potrad(366))
!allocte space for tdew and other variables
allocate(tdew(ndays))
allocate(tmax(ndays))
allocate(tmin(ndays))
allocate(s_tskc(ndays))
!
allocate(optam(21))
optam =(/2.90,3.05,3.21,3.39,3.69,3.82,4.07,4.37,4.72,5.12,5.60,6.18,6.88,7.77,8.90,&
10.39,12.44,15.36,19.79,26.96,30.00/)
!
!----------------------------------------
! estimate srad and humidity with iterative algorithm
! without Tdew input data, an iterative estimation of shortwave radiation and humidity is required
!
! move intent(IN) variables to be able to write on
tmax = s_tmax
tmin = s_tmin
! calculate diurnal temperature range for transmittance calculations
do i=1, ndays
   if (tmax(i) .lt. tmin(i)) tmax(i) = tmin(i)
   dtr(i) = tmax(i) - tmin(i)
end do
!
! smooth dtr array: After Bristow and Campbell, 1984
if (ndays .ge. 30) then !use 30-day antecedent smoothing window
   sm_dtr = pulled_boxcar(dtr=dtr, n=ndays, w=30, w_flag=0)
else !smoothing window width = ndays
   sm_dtr = pulled_boxcar(dtr=dtr, n=ndays, w=ndays, w_flag=0)
end if
! calculate the annual total precip for decision between simple
! and arid-corrected humidity algorithm
sum_prcp = 0.0
do i=1, ndays
   sum_prcp = sum_prcp + s_prcp(i)
end do
ann_prcp = (sum_prcp / ndays) *365.25
if (ann_prcp .eq. 0.0) ann_prcp = 1.0
!
! Generate the effective annual precip, based on a 3-month moving-window.
! Requires some special case handling for the beginning of the record and for short records.
! check if there are at least 90 days in this input file, if not,
! use a simple total scaled to effective annual precip
if (ndays .lt. 90) then
   sum_prcp = 0.0
   do i=1, ndays
      sum_prcp = sum_prcp + s_prcp(i)
   end do
   effann_prcp = (sum_prcp / ndays) *365.35
   ! if the effective annual precip for this period is less than 8 cm,
   ! set the effective annual precip to 8 cm to reflect an arid condition,
   ! while avoiding possible division-by-zero errors and very large ratios (PET/Pann)
   if (effann_prcp .lt. 8.0) effann_prcp = 8.0
   do i=1, ndays
      parray(i) = effann_prcp
   end do
else
   ! Check if the yeardays at beginning and the end of this input file match up.
   ! If so, use parts of the three months at the end of the input file to generate
   ! effective annual precip for the first 3-months.
   ! Otherwise, duplicate the first 90 days of the record
   start_yday = yday(1)
   end_yday = yday(ndays)
   !
   if (.not.(start_yday .eq. 1_i4)) then
      if (end_yday .eq. start_yday-1) then
         isloop = 1
      else
         isloop = 0
      end if
   else
      if ((end_yday .eq. 365) .or. (end_yday .eq. 366)) then
         isloop = 1
      else
         isloop = 0
      end if
   end if
   ! fill the first 90 days of window
   do i=1, 90
      if (isloop .eq. 1_i4) then
         window(i) = s_prcp(ndays-90+i)
      else
         window(i) = s_prcp(i)
      end if
   end do
   ! fill the rest of the window array
   do i=1, ndays
      window(i+90) = s_prcp(i)
   end do
   ! for each day, calculate the effective annual precip from scaled 90-day total
   do i=1, ndays
      sum_prcp = 0.0
      do j=0, 89
         sum_prcp = sum_prcp + window(i+j)
      end do
      sum_prcp = (sum_prcp / 90.0) *365.25
      ! if the effective annual precip for this 90-day period is less than 8 cm,
      ! set the effective annual precip to 8 cm to reflect an arid condition,
      ! while avoiding possible division-by-zero errors and very large ratios (PET/Pann)
      if (sum_prcp .lt. 8.0) then
         parray(i) = 8.0
      else
         parray(i) = sum_prcp
      end if
      !print *, parray(i)
   end do
end if ! end if ndays >= 90
!
!   	 *****************************************
!	 *                                       *
!	 * start of the main radiation algorithm *
!	 *                                       *
!	 *****************************************
!
! before starting the iterative algorithm between humidity and radiation,
! calculate all the variables that don't depend on humidity so they only get done once.
!
! STEP (1) calculate pressure ratio (site/reference) = f(elevation)
t1 = 1.0 - (LR_STD * site_elev) / T_STD
t2 = Gravity_dp / (LR_STD * (R/MA))
pratio = t1 ** t2
!
! STEP (2) correct initial transmittance for elevation
trans1 = TBASE ** pratio
!
! STEP (3) build 366-day array of ttmax0, potential rad, and daylength
! precalculate the transcendentals
lat = site_lat
!
! check for (+/-) 90 degrees latitude, throws off daylength calc
lat = lat * deg2rad_dp
if (lat .gt. 1.5707) lat = 1.5707
if (lat .lt. -1.5707) lat = -1.5707
coslat = cos(lat)
sinlat = sin(lat)
cosslp = cos(site_slp * deg2rad_dp)
sinslp = sin(site_slp * deg2rad_dp)
cosasp = cos(site_asp * deg2rad_dp)
sinasp = sin(site_asp * deg2rad_dp)
!
! cosine of zenith angle for east and west horizons
coszeh = cos(1.570796 - (site_ehoriz * deg2rad_dp))
coszwh = cos(1.570796 - (site_whoriz * deg2rad_dp))
!
! sub-daily time and angular increment information
dt = 600._dp            ! set timestep
dh = dt / SECPERRAD     ! calculate hour-angle step */
!
! begin loop through yeardays
do i=1, 365
   !calculate cos and sin of declination
   decl = MINDECL * cos((i - 1.0 + DAYSOFF) * RADPERDAY)
   cosdecl = cos(decl)
   sindecl = sin(decl)
   ! do some precalculations for beam-slope geometry (bsg)
   bsg1 = -sinslp * sinasp * cosdecl
   bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl
   bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl
   ! calculate daylength as a function of lat and decl
   cosegeom = coslat * cosdecl
   sinegeom = sinlat * sindecl
   coshss = -(sinegeom) / cosegeom
   if (coshss .lt. -1.0) coshss = -1.0  ! 24-hr daylight
   if (coshss .gt.  1.0) coshss =  1.0  ! 0-hr daylight
   hss = acos(coshss)                   ! hour angle at sunset (radians) */
   ! daylength (seconds)
   daylength(i) = 2.0 * hss * SECPERRAD
   ! solar constant as a function of yearday (W/m^2)
   sc = 1368.0 + 45.5 *sin((2.0 *PI_dp *i / 365.25) + 1.7)
   ! extraterrestrial radiation perpendicular to beam, total over the timestep (J)
   dir_beam_topa = sc * dt
   !
   sum_trans = 0.0
   sum_flat_potrad = 0.0
   sum_slope_potrad = 0.0
   !
   ! begin sub-daily hour-angle loop, from -hss to hss
   do h=-hss, hss, dh
   ! precalculate cos and sin of hour angle
      cosh_ = cos(h)
      sinh_ = sin(h)
      ! calculate cosine of solar zenith angle
      cza = cosegeom * cosh_ + sinegeom
      ! calculate cosine of beam-slope angle
      cbsa = sinh_ * bsg1 + cosh_ * bsg2 + bsg3
      ! check if sun is above a flat horizon
      if (cza .gt. 0.0) then
      ! when sun is above the ideal (flat) horizon,
      ! do all the flat-surface calculations to determine daily total
      ! transmittance, and save flat-surface potential radiation for
      ! later calculations of diffuse radiation
      ! potential radiation for this time period, flat surface, top of atmosphere
         dir_flat_topa = dir_beam_topa * cza
      ! determine optical air mass
         am = 1.0 / (cza + 0.0000001)
         if (am .gt. 2.9) then
            ami = int(acos(cza)/deg2rad_dp) - 69 + 1
            if (ami .le. 0_i4) ami = 1
            if (ami .gt. 20_i4) ami = 21
            am = optam(ami)
         end if
         ! correct instantaneous transmittance for this optical air mass
         trans2 = trans1 ** am
         ! instantaneous transmittance is weighted by potential radiation for flat surface at
         ! top of atmosphere to get daily total transmittance
         sum_trans = sum_trans + trans2 * dir_flat_topa
         ! keep track of total potential radiation on a flat surface for ideal horizons
         sum_flat_potrad = sum_flat_potrad + dir_flat_topa
         ! keep track of whether this time step contributes to component 1 (direct on slope)
         ! sun between east and west horizons, and direct on slope.
         ! this period contributes to component 1
         if (((h .lt. 0.0) .and. (cza .gt. coszeh) .and. (cbsa .gt. 0.0)) .or. &
              ((h .ge. 0.0) .and. (cza .gt. coszwh) .and. (cbsa .gt. 0.0))) then
            sum_slope_potrad = sum_slope_potrad + dir_beam_topa * cbsa
         end if
      end if ! end if sun above ideal horizon
   end do ! end of sub-daily hour-angle loop
   !
   ! calculate maximum daily total transmittance and daylight average flux density
   ! for a flat surface and the slope
   if (daylength(i) .gt. 0.0) then
   ! c original: if (daylength[i])
      ttmax0(i) = sum_trans / sum_flat_potrad
      flat_potrad(i) = sum_flat_potrad / daylength(i)
      slope_potrad(i) = sum_slope_potrad / daylength(i)
   else
      ttmax0(i) = 0.0
      flat_potrad(i) = 0.0
      slope_potrad(i) = 0.0
   end if
end do ! end of i=365 days loop
!
! force yearday 366 = yearday 365
ttmax0(366) = ttmax0(365)
flat_potrad(366) = flat_potrad(365)
slope_potrad(366) = slope_potrad(365)
!
! STEP (4)  calculate the sky proportion for diffuse radiation
! uses the product of spherical cap defined by average horizon angle
! and the great-circle truncation of a hemisphere. this factor does not
! vary by yearday.
avg_horizon = (site_ehoriz + site_whoriz) / 2.0
horizon_scalar = 1.0 - sin(avg_horizon * deg2rad_dp)
!
if (site_slp .gt. avg_horizon) then
   slope_excess = site_slp - avg_horizon
else
   slope_excess = 0.0
end if
!
if (2.0 *avg_horizon .gt. 180.0) then
   slope_scalar = 0.0
else
   slope_scalar = 1.0 - (slope_excess / (180.0 - 2.0* avg_horizon))
   if (slope_scalar .lt. 0.0) slope_scalar = 0.0
end if
!
sky_prop = horizon_scalar * slope_scalar
!
! b parameter, and t_fmax not varying with Tdew, so these can be
! calculated once, outside the iteration between radiation and humidity estimates.
! Requires storing t_fmax in an array.
do i=1, ndays
   ! b parameter from 30-day average of DTR
   b = B0 + B1 * exp(-B2 * sm_dtr(i))
   ! proportion of daily maximum transmittance
   t_fmax(i) = 1.0 - 0.9 * exp(-b * dtr(i) **C)
   ! correct for precipitation if this is a rain day
   if (s_prcp(i) .gt. 0.0) t_fmax(i) = t_fmax(i) *RAIN_SCALAR
   ! c original: if (data->prcp[i])
end do
!
! As a first approximation, calculate radiation assuming that Tdew = Tmin
do i=1, ndays
   j = yday(i)
   ! c original: yday = data->yday[i]-1
   tdew(i) = s_tmin(i)
   pva = 610.7 * exp(17.38 * tdew(i) / (239.0 + tdew(i)))
   t_tmax = ttmax0(j) + ABASE * pva
   ! final daily total transmittance
   t_final = t_tmax * t_fmax(i)
   ! estimate fraction of radiation that is diffuse, on an
   ! instantaneous basis, from relationship with daily total
   ! transmittance in Jones (Plants and Microclimate, 1992)
   ! Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
   ! Fig 6.14, p. 122.
   pdif = -1.25 *t_final + 1.25
   if (pdif .gt. 1.0) pdif = 1.0
   if (pdif .lt. 0.0) pdif = 0.0
   !
   ! estimate fraction of radiation that is direct, on an instantaneous basis
   pdir = 1.0 - pdif
   !
   ! the daily total radiation is estimated as the sum of the following two components:
   ! 1. The direct radiation arriving during the part of
   !    the day when there is direct beam on the slope.
   ! 2. The diffuse radiation arriving over the entire daylength
   !    (when sun is above ideal horizon).
   !
   ! component 1
   srad1 = slope_potrad(j) * t_final * pdir
   !
   ! component 2 (diffuse)
   ! includes the effect of surface albedo in raising the diffuse
   ! radiation for obstructed horizons
   srad2 = flat_potrad(j) *t_final *pdif *(sky_prop + DIF_ALB *(1.0 - sky_prop))
   !
   ! snow pack influence on radiation / not sure how to imlement
   if (s_swe(i) .gt. 0.0) then
      ! snow correction in J/m2/day
      sc = (1.32 + 0.096 * s_swe(i)) * 1e6
      ! convert to W/m2 and check for zero daylength
      if (daylength(j) .gt. 0.0) then
         sc = sc / daylength(j)
      else
         sc = 0.0
      end if
      ! set a maximum correction of 100 W/m2
      if (sc .gt. 100.0) sc = 100.0
   else
      sc = 0.0
   end if
   ! save daily radiation and daylength
   s_srad(i) = srad1 + srad2 + sc
   s_srad_dif(i) = srad2
   s_dayl(i) = daylength(j)
end do
!
! estimate annual PET first, to decide which humidity algorithm should be used
! estimate air pressure at site
pa = atm_pres(elev=site_elev)
sum_pet = 0.0
do i=1, ndays
   save_pet(i) = calc_pet(s_srad(i), s_tday(i), pa, s_dayl(i))!, CP, EPS)
   sum_pet = sum_pet + save_pet(i)
end do
ann_pet = (sum_pet/ndays) * 365.25
!
! humidity algorithm decision:
! PET/prcp >= 2.5 -> arid correction
! PET/prcp <  2.5 -> use tdew-tmin, which is already finished
print *, "PET/PRCP = ", ann_pet/ann_prcp
!
if (ann_pet/ann_prcp .ge. 2.5) then
   print *, "Using arid-climate humidity algorithm"
   !
   ! Estimate Tdew using the initial estimate of radiation for PET
   do i=1, ndays
      tmink = s_tmin(i) + 273.15
      pet = save_pet(i)
      !
      ! calculate ratio (PET/effann_prcp) and correct the dewpoint
      ratio = pet / parray(i)
      ratio2 = ratio *ratio
      ratio3 = ratio2 *ratio
      tdewk = tmink *(-0.127 + 1.121 *(1.003 - 1.444*ratio + 12.312*ratio2 - &
           32.766*ratio3) + 0.0006*(dtr(i)))
      tdew(i) = tdewk - 273.15
   end do
   !
   ! Revise estimate of radiation using new Tdew
   do i=1, ndays
      j = yday(i)
      pva = 610.7 * exp(17.38 * tdew(i) / (239.0 + tdew(i)))
      t_tmax = ttmax0(j) + ABASE * pva
      !
      ! final daily total transmittance
      t_final = t_tmax * t_fmax(i)
      !
      ! estimate fraction of radiation that is diffuse, on an
      ! instantaneous basis, from relationship with daily total
      ! transmittance in Jones (Plants and Microclimate, 1992)
      ! Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
      ! Fig 6.14, p. 122.
      pdif = -1.25 *t_final + 1.25
      if (pdif > 1.0) pdif = 1.0
      if (pdif < 0.0) pdif = 0.0
      !
      ! estimate fraction of radiation that is direct, on an instantaneous basis
      pdir = 1.0 - pdif
      !
      ! the daily total radiation is estimated as the sum of the following two components:
      ! 1. The direct radiation arriving during the part of
      !    the day when there is direct beam on the slope.
      ! 2. The diffuse radiation arriving over the entire daylength
      !    (when sun is above ideal horizon).
      !
      ! component 1
      srad1 = slope_potrad(j) * t_final * pdir
      ! component 2 (diffuse)
      ! includes the effect of surface albedo in raising the diffuse
      ! radiation for obstructed horizons
      srad2 = flat_potrad(j) *t_final *pdif *(sky_prop + DIF_ALB *(1.0-sky_prop))
      ! snow pack influence on radiation / not sure how to implement
      if (s_swe(i) .gt. 0.0) then
         ! snow correction in J/m2/day
         sc = (1.32 + 0.096 * s_swe(i)) * 1e6
         ! convert to W/m2 and check for zero daylength
         if (daylength(j) .gt. 0.0) then
            sc = sc / daylength(j)
         else
            sc = 0.0
         end if
         ! set a maximum correction of 100 W/m2
         if (sc .gt. 100.0) sc = 100.0
      else
         sc = 0.0
      end if
      ! save daily radiation and daylength
      s_srad(i) = srad1 + srad2 + sc
      s_srad_dif(i) = srad2
      !daylength(i) = daylength(j)
   end do
   !
   ! Revise estimate of Tdew using new radiation
   do i=1, ndays
      tmink = s_tmin(i) + 273.15
      pet = save_pet(i)
      !
      ! calculate ratio (PET/effann_prcp) and correct the dewpoint
      ratio = pet / parray(i)
      ratio2 = ratio *ratio
      ratio3 = ratio2 *ratio
      tdewk = tmink *(-0.127 + 1.121 *(1.003 - 1.444*ratio + 12.312*ratio2 - &
           32.766*ratio3) + 0.0006*(dtr(i)))
      tdew(i) = tdewk - 273.15
   end do
else       ! end of arid-correction humidity and radiation estimation
   print *, "Using Tdew=Tmin humidity algorithm"
end if
!
! now calculate vapor pressure from tdew
do i=1, ndays
   pva = 610.7 * exp(17.38 * tdew(i) / (239.0 + tdew(i)))
   select case (outhum)
      case (0)
        s_hum(i) = pva
      case (1)
        ! output humidity as vapor pressure deficit (Pa)
        ! calculate saturated VP at tday
        pvs = 610.7 * exp(17.38 * s_tday(i) / (239.0 + s_tday(i)))
        vpd = pvs - pva
        if (vpd .lt. 0.0) vpd = 0.0
        s_hum(i) = vpd
    end select
end do
!
!*/begin longwave
do i=1, ndays
   j = yday(i)
   select case (lwrad)
      case (0)
         ! no calc of longwave radiation
         s_lrad(i) = -9999.
      case (1)
         ! Longwave correction factor, used to correct estimated incoming longwave radiation
         ! (use 1, unless measured longwave available for calibration)
         ! LWAVE_COR = parameterset(14)
         ! cloud cover fraction of VIC implementation mtCLIM
         ! /* See Bras, R. F. , "Hydrology, an introduction to hydrologic science", Addison-Wesley, 1990, p. 42-45 */
         ! /* start vic_change */
         pva = 610.7 * exp(17.38 * tdew(i) / (239.0 + tdew(i)))
         s_tskc(i) = sqrt((1. - t_fmax(i))/ 0.65)
         emissivity = 0.740 + 0.0049 * pva / 100.
         cloudfactor = 1.0 + 0.17 * s_tskc(i) * s_tskc(i)
         s_lrad(i) = emissivity * cloudfactor * sigma_dp *(s_tday(i) + 273.15_dp)**4.0_dp / LWAVE_COR
         ! /* end vic_change */
      case (2)
         pva = 610.7 * exp(17.38 * tdew(i) / (239.0 + tdew(i)))
         ! Bohn et al. 2013
         ! clear sky shortwave radiation with t_tmax
         ! component 1
         srad1 = slope_potrad(j) * t_tmax * pdir
         ! component 2 (diffuse)
         ! includes the effect of surface albedo in raising the diffuse
         ! radiation for obstructed horizons
         srad2 = flat_potrad(j) *t_tmax *pdif *(sky_prop + DIF_ALB *(1.0-sky_prop))
         ! source emissivity for full sky (Deardorff) is 1
         ! source emissivity E_a for clear sky (Prata)
         ! convert actual vapor pressure (vpa) from Pa to hPa
         ! convert mean air temperature from degC to Kelvin
         xi = 46.5_dp * pva/100.0_dp / (s_tday(i) + 273.15_dp)
         E_a = 1.0_dp - (1.0_dp+xi)* exp(-sqrt(1.2_dp+3.0_dp*xi))
         ! incoming longwave radiation
         ! fraction of clear sky: s_srad(i)/slope_potrad(i) - actual shortwave radiation / potential shortwave radiation
         s_lrad(i) = ((1.0_dp - s_srad(i)/(srad1 + srad2 + sc)) + s_srad(i)/(srad1 + srad2 + sc)* E_a) &
              *sigma_dp *(s_tday(i) + 273.15_dp)**4.0_dp
    end select

    select case (netlwrad)
      case (0)
        ! no calc of net longwave radiation
        s_lradnet(i) = -9999.
      case (1)
        ! outgoing net longwave radiation
        ! Longwave correction factor, used to correct estimated incoming longwave radiation
        ! (use 1, unless measured longwave available for calibration)
        ! LWAVE_COR = parameterset(14)
        ! cloud cover fraction of VIC implementation mtCLIM
        ! /* See Bras, R. F. , "Hydrology, an introduction to hydrologic science", Addison-Wesley, 1990, p. 42-45 */
        pva = 610.7 * exp(17.38 * tdew(i) / (239.0 + tdew(i)))
        s_tskc(i) = sqrt((1. - t_fmax(i))/ 0.65)
        emissivity = 0.740 + 0.0049 * pva / 100.
        cloudfactor = 1.0 + 0.17 * s_tskc(i) * s_tskc(i)
        s_lradnet(i) = (-1.) * (emissivity * cloudfactor / LWAVE_COR - 1 )* sigma_dp* (s_tday(i)+ 273.15_dp)**4
      case(2)
        pva = 610.7 * exp(17.38 * tdew(i) / (239.0 + tdew(i)))
        ! clear sky shortwave radiation with t_tmax
        ! component 1
        srad1 = slope_potrad(j) * t_tmax * pdir
        ! component 2 (diffuse)
        ! includes the effect of surface albedo in raising the diffuse
        ! radiation for obstructed horizons
        srad2 = flat_potrad(j) *t_tmax *pdif *(sky_prop + DIF_ALB *(1.0-sky_prop))
        ! outgoing net longwave radiation (Allen et al. 1998, FAO)
        s_lradnet(i) = sigma_dp* (s_tday(i)+ 273.15_dp)**4 * &
             (0.34_dp-0.14_dp* sqrt(pva/1000.0_dp))* (1.35_dp* s_srad(i)/(srad1 + srad2 + sc) - 0.35_dp )
    end select

end do
!*/end longwave

end subroutine calc_srad_humidity_iterative


  ! ------------------------------------------------------------------------------
  !    NAME
  !        pulled_boxcar
  !
  !    PURPOSE
  !        calculates a moving average of antecedent values in an array,
  !        using either a ramped (w_flag=1) or a flat (w_flag=0) weighting

  !    CALLING SEQUENCE
  !        sm_dtr = pulled_boxcar(dtr, n, w, w_flag)

  !    INTENT(IN)
  !        character(len=*) :: dtr - diurnal temperature range

  !    INTENT(OUT)
  !        integer(i4)      :: n - number of yearday values

  !    INTENT(IN)
  !        integer(i4)      :: w - weighting length

  !    INTENT(IN)
  !        integer(i4)      :: w_flag - switch between ramped (w_flag=1) or
  !                                     flat weighting (w_flag=0)

  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        none

  !    LITERATURE
  !        Bristow and Campbell, 1984

  !    HISTORY
  !        Written,  Johannes Brenner, Sept 2016

  function pulled_boxcar (dtr, n, w, w_flag)
    !
    real(dp),     dimension(:), allocatable, intent(in)   :: dtr
    real(dp),     dimension(:), allocatable               :: pulled_boxcar
    !
    integer(i4)                            , intent(in)   :: n
    integer(i4)                            , intent(in)   :: w
    integer(i4)                            , intent(in)   :: w_flag
    !
    integer(i4)                                           :: ok =1
    integer(i4)                                           :: i, j
    real(dp),     dimension(:), allocatable               :: wt
    real(dp)                                              :: total, sum_wt
    !
    ! allocate
    allocate(pulled_boxcar(n))
    allocate(wt(w))
    !
    if (w .gt. n) then
        print *, "Boxcar longer than array...\n"
        print *, "Resize boxcar and try again\n"
        ok=0
    end if
    !
    if (ok .eq. 1_i4) then
    ! if w_flag != 0, use linear ramp to weight tails, otherwise use constant weight
      sum_wt = 0.0
      if (w_flag .eq. 1_i4) then
         do i=1, w
            wt(i) = i
            sum_wt = sum_wt + wt(i)
         end do
      else
         do i=1, w
            wt(i) = 1.0
            sum_wt = sum_wt + wt(i)
         end do
      end if
    !
    ! fill the output array, starting with the point where a full boxcar can be calculated
      do i=w, n
         total = 0.0
         do j=1, w
            total = total + dtr(i-w+j) * wt(j)
         end do
         pulled_boxcar(i) = total / sum_wt
      end do
    !
    ! fill the first w elements of the output array with the value from the first full boxcar
      do i=1, w-1
         pulled_boxcar(i) = pulled_boxcar(w)
      end do
    !
    end if
  end function pulled_boxcar

  ! ------------------------------------------------------------------------------
  !    NAME
  !        atm_pres
  !
  !    PURPOSE
  !        calculates the atmospheric pressure as a function of elevation

  !    CALLING SEQUENCE
  !        sm_dtr = atm_pres(elev)

  !    INTENT(IN)
  !        real(dp) :: elev - elevation of site (m)
  !

  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        none

  !    LITERATURE
  !        Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics, 2nd
  !		Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
  !		(p. 168)
  !    HISTORY
  !        Written,  Johannes Brenner, Sept 2016

function atm_pres(elev)
   !
   ! daily atmospheric pressure (Pa) as a fu nction of elevation (m)
   real(dp), intent(in)   :: elev
   real(dp)               :: atm_pres
   real(dp)               :: t1, t2
   !
   ! (-K m-1) standard temperature lapse rate
   !real(dp)               :: LR_STD = 0.0065
   ! (K) standard temp at 0.0 m elevation
   !real(dp)               :: T_STD = 288.15
   ! (m s-2) standard gravitational accel. - Gravity_dp
   !real(dp)               :: G_STD = 9.80665
   ! (Pa) standard pressure at 0.0 m elevation
   !real(dp)               :: P_STD = 101325.0 - P0_dp
   ! (m3 Pa mol-1 K-1) gas law constant
   !real(dp)               :: R = 8.3143
   ! (kg mol-1) molecular weight of air
   !real(dp)               :: MA = 28.9644e-3

   !
   t1 = 1.0 - (LR_STD * elev)/T_STD
   t2 = Gravity_dp / (LR_STD * (R / MA))
   atm_pres = P0_dp * (t1 **t2)
   !
end function atm_pres

  ! ------------------------------------------------------------------------------
  !    NAME
  !        snowpack
  !
  !    PURPOSE
  !        estimates the accumulation and melt of snow for radiation algorithm corrections

  !    CALLING SEQUENCE
  !        s_swe = snowpack(tmin, prcp, yday, ndays, SNOW_TCRIT, SNOW_TRATE)

  !    INTENT(IN)
  !        real(dp), dimension(:), allocatable :: tmin - minimum day air temperature of site (degC)
  !
  !        real(dp), dimension(:), allocatable :: prcp - daily precipitation sum of site (mm)
  !
  !        integer(i4) :: yday - yearday values
  !
  !        integer(i4) :: ndays - yearday values
  !
  !        integer(i4) :: yday - number of days of data
  !
  !        integer(i4) :: SNOW_TCRIT - critical temperature for snowmelt (degC)
  !
  !        integer(i4) :: SNOW_TRATE - snowmelt rate (cm/degC/day)
  !
  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        none

  !    LITERATURE
  !
  !    HISTORY
  !        Written,  Johannes Brenner, Sept 2016

function snowpack(tmin, prcp, yday, ndays, SNOW_TCRIT, SNOW_TRATE)
   !
   real(dp),    dimension(:), allocatable, intent(in) :: tmin
   real(dp),    dimension(:), allocatable, intent(in) :: prcp
   integer(i4), dimension(:), allocatable, intent(in) :: yday
   real(dp)                              , intent(in) :: SNOW_TCRIT
   real(dp)                              , intent(in) :: SNOW_TRATE
   real(dp), dimension(:), allocatable                :: snowpack
   !
   real(dp)                    :: snowpack_,newsnow,snowmelt,sum_
   !
   integer(i4)                 :: i, ndays, counter
   integer(i4)                 :: start_yday, prev_yday
   !
   allocate(snowpack(ndays))
   !
   ! first pass to initialize SWE array
   snowpack_ = 0.0
   !
   do i=1, ndays
      newsnow  = 0.0
      snowmelt = 0.0
      if (tmin(i) .le. SNOW_TCRIT) then
         newsnow = prcp(i)
      else
         snowmelt = SNOW_TRATE * (tmin(i) - SNOW_TCRIT)
      end if
      snowpack_ = snowpack_ + newsnow - snowmelt
      if (snowpack_ .lt. 0.0) snowpack_ = 0.0
      snowpack(i) = snowpack_
      !print *, snowpack_
   end do
   ! use the first pass to set the initial snowpack conditions for the
   ! first day of data
   start_yday = yday(1)
   if (start_yday .eq. 1_i4) then
      prev_yday = 365
   else
      prev_yday = start_yday - 1
   end if
   counter = 0
   sum_ = 0.0
   do i=2, ndays
      if ((yday(i) .eq. start_yday) .or. (yday(i) .eq. prev_yday)) then
         counter = counter + 1
         sum_ = sum_ + snowpack(i)
      end if
   end do
   ! Proceed with correction if there are valid days to reinitialize
   ! the snowpack estiamtes. Otherwise use the first-pass estimate.
   if (counter .gt. 0) then
      snowpack_ = sum_ / real(counter)
      do i=1, ndays
         newsnow  = 0.0
         snowmelt = 0.0
         if (tmin(i) .le. SNOW_TCRIT) then
            newsnow = prcp(i)
         else
            snowmelt = SNOW_TRATE * (tmin(i) - SNOW_TCRIT)
         end if
         snowpack_ = snowpack_ + newsnow - snowmelt
         !print *, snowpack_
      if (snowpack_ .lt. 0.0) snowpack_ = 0.0
      snowpack(i) = snowpack_
      end do
   end if
end function snowpack

  ! ------------------------------------------------------------------------------
  !    NAME
  !        calc_pet
  !
  !    PURPOSE
  !        calculates the potential evapotranspiration for aridity
  !        corrections in calc_vpd(), according to Kimball et al., 1997

  !    CALLING SEQUENCE
  !        pet = calc_pet(rad, ta, pa, dayl)

  !    INTENT(IN)
  !        real(dp)    :: rad -  daylight average incident shortwave radiation (W/m2)
  !
  !        real(dp)    :: ta - daylight average air temperature (deg C)
  !
  !        real(dp)    :: pa - air pressure (Pa)
  !
  !        real(dp)    :: dayl - daylength
  !
  !        real(dp)    :: CP - specific heat of air (J/kg K)
  !
  !        real(dp)    :: EPS - unitless ratio of molec weights
  !
  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        none

  !    LITERATURE
  !        Abbott, P.F., and R.C. Tabony, 1985. The estimation of humidity parameters.
  !                            Meteorol. Mag., 114:49-56.
  !        Kimball, J.S., S.W. Running, and R. Nemani, 1997. An improved method for estimating
  !                              surface humidity from daily minimum temperature.
  !                              Agricultural and Forest Meteorology, 85:87-98.
  !

  !    HISTORY
  !        Written,  Johannes Brenner, Sept 2016

! JBJBJBJB
function calc_pet(rad, ta, pa, dayl)!, CP=cp0_dp, EPS)
   !
   real(dp), intent(in) :: rad
   real(dp), intent(in) :: ta
   real(dp), intent(in) :: pa
   real(dp), intent(in) :: dayl
   real(dp)             :: calc_pet ! (kg m-2 day-1) potential evapotranspiration
   !
   real(dp)             :: rnet   ! (W m-2) absorbed shortwave radiation avail. for ET
   real(dp)             :: lhvap  ! (J kg-1) latent heat of vaporization of water
   real(dp)             :: gamma_ ! (Pa K-1) psychrometer parameter
   real(dp)             :: dt = 0.2  ! offset for saturation vapor pressure calculation
   real(dp)             :: t1, t2    ! (deg C) air temperatures
   real(dp)             :: pvs1, pvs2  ! (Pa)   saturated vapor pressures
   real(dp)             :: s      ! (Pa K-1) slope of saturated vapor pressure curve
   !
   ! calculate absorbed radiation, assuming albedo = 0.2  and ground
   ! heat flux = 10% of absorbed radiation during daylight
   rnet = rad * 0.72
   !
   ! calculate latent heat of vaporization as a function of ta
   lhvap = 2.5023e6 - 2430.54 * ta
   !
   ! calculate the psychrometer parameter: gamma = (cp pa)/(lhvap epsilon) where:
   ! cp       (J/kg K)   specific heat of air
   ! epsilon  (unitless) ratio of molecular weights of water and air
   gamma_ = cp0_dp * pa / (lhvap * EPS)
   !
   ! estimate the slope of the saturation vapor pressure curve at ta
   ! temperature offsets for slope estimate
   t1 = ta+dt
   t2 = ta-dt
   !
   ! calculate saturation vapor pressures at t1 and t2, using formula from
   ! Abbott, P.F., and R.C. Tabony, 1985. The estimation of humidity parameters.
   ! Meteorol. Mag., 114:49-56.
   pvs1 = 610.7 * exp(17.38 * t1 / (239.0 + t1))
   pvs2 = 610.7 * exp(17.38 * t2 / (239.0 + t2))
   !
   ! calculate slope of pvs vs. T curve near ta
   s = (pvs1-pvs2) / (t1-t2)
   !
   ! calculate PET using Priestly-Taylor approximation, with coefficient
   ! set at 1.26. Units of result are kg/m^2/day, equivalent to mm water/day
   calc_pet = (1.26 * (s / (s+gamma_)) * rnet * dayl) / lhvap
   !
   ! return a value in centimeters/day, because this value is used in a ratio
   ! to annual total precip, and precip units are centimeters
   calc_pet = calc_pet / 10
end function calc_pet

  ! ------------------------------------------------------------------------------
  !    NAME
  !        calc_tair
  !
  !    PURPOSE
  !        calculates daily min, max, mean air temperature at site

  !    CALLING SEQUENCE
  !        call calc_tair(tmin, tmax, ndays, site_elev, base_elev, TDAYCOEF, tmin_lr, tmax_lr, s_tmax, s_tmin, s_tday)
  !
  !    INTENT(IN)
  !        real(dp), dimension(:), allocatable, intent(in) :: tmin - base minimum day air temperature (degC)
  !
  !        real(dp), dimension(:), allocatable, intent(in) :: tmax - base maximum day air temperature (degC)
  !
  !        integer(i4), intent(in) :: ndays - yearday values
  !
  !        real(dp)   , intent(in) :: site_elev - Site elevation, meters
  !
  !        real(dp)   , intent(in) :: base_elev - Base elevation, meters
  !
  !        real(dp)   , intent(in) :: TDAYCOEF - daylight air temperature coefficient (dim)
  !
  !        real(dp)   , intent(in) :: tmin_lr - Minimum temperature lapse rate (deg C/ km)
  !
  !        real(dp)   , intent(in) :: tmax_lr - Maximum temperature lapse rate (deg C/ km)
  !
  !        real(dp), dimension(:), allocatable, intent(out) :: s_tmin - minimum day air temperature of site (degC)
  !
  !        real(dp), dimension(:), allocatable, intent(out) :: s_tmax - maximum day air temperature of site (degC)
  !
  !        real(dp), dimension(:), allocatable, intent(out) :: s_tday - mean day air temperature of site (degC)
  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        none

  !    LITERATURE
  !

  !    HISTORY
  !        Written,  Johannes Brenner, Sept 2016

subroutine calc_tair(tmin, tmax, ndays, site_elev, base_elev, TDAYCOEF, tmin_lr, tmax_lr,&
     s_tmax, s_tmin, s_tday)
  !
  real(dp),    dimension(:), allocatable, intent(in)   :: tmin
  real(dp),    dimension(:), allocatable, intent(in)   :: tmax
  !
  real(dp),                               intent(in)   :: site_elev
  real(dp),                               intent(in)   :: base_elev
  real(dp),                               intent(in)   :: TDAYCOEF
  real(dp),                               intent(in)   :: tmin_lr
  real(dp),                               intent(in)   :: tmax_lr
  integer(i4),                            intent(in)   :: ndays
  !
  real(dp),    dimension(:), allocatable, intent(out)  :: s_tmin
  real(dp),    dimension(:), allocatable, intent(out)  :: s_tmax
  real(dp),    dimension(:), allocatable, intent(out)  :: s_tday
  !
  real(dp)                                             :: dz
  real(dp)                                             :: tmean
  integer(i4)                                          :: i
  !
  ! allocate
  allocate(s_tday(ndays))
  allocate(s_tmax(ndays))
  allocate(s_tmin(ndays))
  ! calculate elevation difference in kilometers
  dz = (site_elev - base_elev) / 1000.0
  ! apply lapse rate corrections to tmax and tmin
  ! Since tmax lapse rate usually has a larger absolute value than tmi
  ! lapse rate, it is possible at high elevation sites for these corrections
  ! to result in tmin > tmax. Check for that occurrence and force
  ! tmin = corrected tmax - 0.5 deg C.
  do i=1, ndays
    ! lapse rate corrections
    s_tmax(i) = tmax(i) + (dz * tmax_lr)
    s_tmin(i) = tmin(i) + (dz * tmin_lr)
    ! derived temperatures
    tmean = (s_tmax(i) + s_tmin(i)) / 2.0
    s_tday(i) = ((s_tmax(i) - tmean) *TDAYCOEF) + tmean
  end do
  !
end subroutine calc_tair

  ! ------------------------------------------------------------------------------
  !    NAME
  !        calc_prcp
  !
  !    PURPOSE
  !        calculates daily total precipitation at site

  !    CALLING SEQUENCE
  !        s_prcp = calc_prcp(prcp, ndays, s_prcp)
  !
  !    INTENT(IN)
  !        real(dp), dimension(:), allocatable, intent(in) :: prcp - base daily total precipitation (mm)
  !
  !        integer(i4), intent(in) :: ndays - yearday values
  !
  !        real(dp)   , intent(in) :: site_isoh - Site annual precip isohyet, cm
  !
  !        real(dp)   , intent(in) :: base_isoh - Base annual precip isohyet, cm
  !

  !    RESTRICTIONS
  !

  !    EXAMPLE
  !        none

  !    LITERATURE
  !

  !    HISTORY
  !        Written,  Johannes Brenner, Sept 2016

function calc_prcp(prcp, ndays, site_isoh, base_isoh)
  !
  real(dp),    dimension(:), allocatable, intent(in)   :: prcp
  !
  integer(i4),                            intent(in)   :: ndays
  real(dp),                               intent(in)   :: site_isoh
  real(dp),                               intent(in)   :: base_isoh
  !
  real(dp),    dimension(:), allocatable               :: calc_prcp
  !
  real(dp)                                             :: ratio
  integer(i4)                                          :: i
  !
  allocate(calc_prcp(ndays))
  !
  ratio = site_isoh / base_isoh
  !
  do i=1, ndays
     calc_prcp(i) = prcp(i) *ratio
  end do
end function calc_prcp

subroutine eval_mtclim(parameterset, yday, tmin, tmax, prcp, ndays, &
     outhum, lwrad, netlwrad, &
     s_tmax_out, s_tmin_out, s_tday_out, s_prcp_out, &
     s_hum_out, s_srad_out, s_srad_dif_out, s_lrad_out, s_lradnet_out, s_dayl_out, &
     site_elev, base_elev, tmin_lr, tmax_lr, site_isoh, base_isoh, &
     site_lat, site_asp, site_slp, site_ehoriz, site_whoriz)
  !
  real(dp),     dimension(14), intent(in)               :: parameterset
  real(dp),     dimension(:), allocatable, intent(in)   :: tmin, tmax, prcp
  integer(i4),  dimension(:), allocatable, intent(in)   :: yday
  integer(i4),                             intent(in)   :: ndays
  !
  !select cases
  integer(i4),                             intent(in)   :: outhum, lwrad, netlwrad
  !
  real(dp),     dimension(:), allocatable, intent(out) :: s_tmax_out  !/* array of site tmax values */
  real(dp),     dimension(:), allocatable, intent(out) :: s_tmin_out  !/* array of site tmin values */
  real(dp),     dimension(:), allocatable, intent(out) :: s_tday_out  !/* array of site daylight temperature values */
  real(dp),     dimension(:), allocatable, intent(out) :: s_prcp_out  !/* array of site prcp values */
  real(dp),     dimension(:), allocatable, intent(out) :: s_hum_out   !/* array of site humidity values (VPD or VP, Pa) */
  real(dp),     dimension(:), allocatable, intent(out) :: s_srad_out  !/* array of site shortwave radiation values */
  real(dp),     dimension(:), allocatable, intent(out) :: s_srad_dif_out  !/* array of site shortwave diffuse radiation values */
  real(dp),     dimension(:), allocatable, intent(out) :: s_lrad_out  !/* array of site longwave radiation values */
  real(dp),     dimension(:), allocatable, intent(out) :: s_lradnet_out  !/* array of site net longwave radiation values */
  real(dp),     dimension(:), allocatable, intent(out) :: s_dayl_out  !/* array of site daylength values */

  real(dp),     dimension(:), allocatable :: s_swe   !/* array of site snow water equivalent values (cm) */
  real(dp)     :: site_elev, base_elev, tmin_lr, tmax_lr, site_isoh, base_isoh, &
       site_lat, site_asp, site_slp, site_ehoriz, site_whoriz
  !
  ! allocate space in the data arrays for input and output data
  !
  allocate(s_swe(ndays))
  !
  allocate(s_tmin_out(ndays))
  allocate(s_tmax_out(ndays))
  allocate(s_tday_out(ndays))
  !
  allocate(s_hum_out(ndays))
  allocate(s_srad_out(ndays))
  allocate(s_srad_dif_out(ndays))
  allocate(s_lrad_out(ndays))
  allocate(s_lradnet_out(ndays))
  allocate(s_dayl_out(ndays))
  !
  !----------------------------------------
  ! inital progress indicator
  print *, "Starting MTCLIM version 4.3\n"
  !
  !----------------------------------------
  ! estimate daily air temperatures
  call calc_tair(tmin, tmax, ndays, site_elev, base_elev, parameterset(13), &
                 tmin_lr, tmax_lr, s_tmax_out, s_tmin_out, s_tday_out)
  print *, "Completed calc_tair()"
  !
  !----------------------------------------
  ! estimate daily precipitation
  s_prcp_out = calc_prcp(prcp, ndays, site_isoh, base_isoh)
  print *, "Completed calc_prcp()"
  !----------------------------------------
  ! estimate daily snowpack
  s_swe = snowpack(s_tmin_out, s_prcp_out, yday, ndays, &
                   parameterset(11), parameterset(12))
  print *, "Completed snowpack()"
  !
  !----------------------------------------
  ! estimate srad and humidity with iterative algorithm
  ! without Tdew input data, an iterative estimation of shortwave radiation and humidity is required
  call calc_srad_humidity_iterative( &
     yday, s_tmax_out, s_tmin_out, s_tday_out, s_prcp_out, s_swe, &
     s_hum_out, s_srad_out, s_srad_dif_out, s_lrad_out, s_lradnet_out, s_dayl_out, &
     ndays, outhum, lwrad, netlwrad, &
     site_lat, site_elev, site_asp, site_slp, site_ehoriz, site_whoriz, &
     parameterset(1), parameterset(2), parameterset(3), &
     parameterset(4), parameterset(5), &
     parameterset(6), parameterset(7), parameterset(8), parameterset(14))
  print *, "Completed radiation algorithm"

end subroutine eval_mtclim

end module mo_mtclim
