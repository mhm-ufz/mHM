!> \file mo_neutrons.f90
!> \brief   \copybrief mo_neutrons
!> \details \copydetails mo_neutrons

!> \brief Models to predict neutron intensities above soils
!> \details The number of neutrons above the ground is directly related to
!! the number soil water content in the ground, air, vegetation and/or snow.
!! This module forward-models neutron abundance as a state variable for each cell.
!> \authors Martin Schroen
!> \date Mar 2015
!> \warning THIS MODULE IS WORK IN PROGRESS, DO NOT USE FOR RESEARCH.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_neutrons

  ! TODO make it faster with pre-calculated horizons and variables
  ! TODO use global parameters as linear model

  USE mo_kind, ONLY: i4, dp
  IMPLICIT NONE

  ! Neutron forward model using particle transport physics, see Shuttleworth et al. 2013
  PUBLIC :: COSMIC

  ! inverse \theta(N) relation based on Desilets et al. 2010
  PUBLIC :: DesiletsN0

  ! integration tabular for approximating the neutron flux integral
  PUBLIC :: TabularIntegralAFast

  PRIVATE

CONTAINS

  ! -----------------------------------------------------------------------------------
  !     NAME
  !         DesiletsN0
  !
  !     PURPOSE
  !>        \brief Calculate neutrons from soil moisture for effective soil layer
  !>        \details Using the N0-relation derived by Desilets, neutron
  !>        counts above the ground (one value per cell in mHM) can be
  !>        derived by a semi-empirical, semi-physical relation.
  !>        The result depends on N0, the neutron counts for 0% soil mositure.
  !>        This variable is site-specific and is a global parameter in mHM.
  !         ------------------------------------------------------------------
  !         N0 formula based on Desilets et al. 2010
  !         ------------------------------------------------------------------
  !
  !     CALLING SEQUENCE
  !         call DesiletsN0( Moisture(cells,layers), Depths(layers), &
  !                          N0-parameter, output(cells) )
  !
  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)   :: SoilMoisture" Soil Moisture
  !>        \param[in] "real(dp), dimension(:)   :: Horizon_depth" Horizon depths
  !>        \param[in] "real(dp), dimension(:)   :: Bd"       Bulk density
  !>        \param[in] "real(dp), dimension(:)   :: LatWater" Lattice water
  !>        \param[in] "real(dp)                 :: N0"       dry neutron counts
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(size(SoilMoisture,1)) :: neutrons" Neutron counts
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         Horizons(1) must not be zero.
  !
  !     EXAMPLE
  !         N0=1500cph, SoilMoisture(1,1)=70mm, Horizons(1)=200mm
  !         1500*(0.372+0.0808/ (70mm/200mm + 0.115))
  !         DesiletsN0 = 819cph
  !
  !     LITERATURE
  !         Desilets, D., M. Zreda, and T. P. A. Ferre (2010),
  !         Nature's neutron probe: Land surface hydrology at an elusive scale
  !         with cosmic rays, WRR, 46, W11505, doi:10.1029/2009WR008726.
  !
  !     HISTORY
  !>        \author Martin Schroen
  !>        \date Mar 2015
  !>        Modified Rohini Kumar         Oct 2021 - Vertical weighting approach
  !>                                      for the Neutron count module to mHM - develop branch

  subroutine DesiletsN0(SoilMoisture, Horizon_depth, Bd, latWater, N0, neutrons)

    use mo_mhm_constants, only: Desilets_a0, Desilets_a1, Desilets_a2
    use mo_moment,        only: average
    implicit none

    real(dp), dimension(:),          intent(in)    :: SoilMoisture
    real(dp), dimension(:),          intent(in)    :: Horizon_depth
    real(dp), dimension(:),          intent(in)    :: Bd
    real(dp), dimension(:),          intent(in)    :: latWater
    real(dp),                        intent(in)    :: N0           ! from global parameters
    real(dp),                        intent(inout) :: neutrons
    ! local variables
    integer(i4)                                    :: nLayers, LL, nn, nIntervals
    real(dp), dimension(:), allocatable            :: Layer_min, Layer_max, Layer_depth
    real(dp)                                       :: average_swc, average_bd
    real(dp)                                       :: D_86_in_cm,  D_86_in_mm
    real(dp), dimension(:), allocatable            :: cummulative_Layer_weight
    real(dp)                                       :: depth, weight_10mm_spacing, grav_swc

    ! get the # of soil hoizons
    nLayers = size(SoilMoisture)
    allocate(  Layer_min(nLayers)  )
    allocate(  Layer_max(nLayers)  )
    allocate( Layer_depth(nLayers) )

    ! assign layer-1
    Layer_min(1) = 0.0_dp
    Layer_max(1) = Horizon_depth(1)
    ! rest layers
    do LL = 2, nLayers
       Layer_min(LL) = Layer_max(LL-1)
       Layer_max(LL) = Horizon_depth(LL)
    end do

    ! estimate layer depth [mm]
    Layer_depth(:) = Layer_max(:) - Layer_min(:)

    ! average soil water content (volumetric ones) and Bulk density
    average_swc = average( SoilMoisture(:)/Layer_depth(:) )
    average_bd  = average( Bd(:) )

    ! estimate D86 [in cm and mm ]
    D_86_in_cm = ( 1.0_dp/average_bd ) * &
         ( 8.321_dp + 0.14249_dp * (0.96655_dp + exp(-0.01_dp)) * (20.0_dp + average_swc)/(0.0429_dp + average_swc) )
    D_86_in_mm = D_86_in_cm * 10.0_dp  !# convert cm to mm


    ! initalise cummulative layer specific weight
    allocate( cummulative_Layer_weight(nLayers)  )
    cummulative_Layer_weight(:) = 0.0_dp

    !! only once to calculate weight at equal space at every 10 mm spacing
    nIntervals = nint( maxval(Horizon_depth(:))/10.0_dp )

    ! calculate 10 mm spacing and on fly the cummulative weights
    depth = 0.0_dp
    do nn = 1, nIntervals
       weight_10mm_spacing = exp(-2.0_dp * depth / D_86_in_mm )

       ! estimate the layer specific cummulative weights
       do LL = 1, nLayers
          if(  (depth .GE. Layer_min(LL))  .AND. (depth .LT. Layer_max(LL))  ) then
             cummulative_Layer_weight(LL) = cummulative_Layer_weight(LL) + weight_10mm_spacing
          end if
       end do
       ! update depth
       depth = depth + 10.0_dp
    end do


    ! estimate weightage SWC
    ! here average_swc is gravemtric swc = vol.swc/Bulk density
    average_swc = 0.0_dp
    do LL = 1, nLayers
       !! >> add here if organic water is coming !! >> SoilMoisture(LL) + latWater(LL) + organic_water(LL)...
       grav_swc = ( SoilMoisture(LL) + latWater(LL) )/Layer_depth(LL)/Bd(LL)
       average_swc = average_swc + ( grav_swc * cummulative_Layer_weight(LL) )
    end do
    average_swc = average_swc / sum( cummulative_Layer_weighT(:) )

    ! calculate neutron count based on depth weighted SM of *D86*
    neutrons = N0 * ( Desilets_a1 + Desilets_a0 / (average_swc + Desilets_a2) )

    !! deallocate variables
    deallocate(Layer_min, Layer_max, Layer_depth, cummulative_Layer_weight)

  end subroutine DesiletsN0
  ! -----------------------------------------------------------------------------------
  !     NAME
  !         COSMIC
  !
  !     PURPOSE
  !>        \brief Calculate neutrons from soil moisture in all layers.
  !>        \details Neutron counts above the ground (one value per cell in mHM)
  !>        can be derived by a simplified physical neutron transport simulation.
  !>        Fast cosmic-Ray neutrons are generated in the soil and attenuated
  !>        differently in water and soil. The remaining neutrons that reached
  !>        the surface relate to the profile of soil water content below.
  !>        Variables like N, alpha and L3 are site-specific and need to be calibrated.
  !         ------------------------------------------------------------------
  !         COSMIC model based on Shuttleworth et al. 2013
  !         ------------------------------------------------------------------
  !
  !     CALLING SEQUENCE
  !         call COSMIC( Moisture(cells,layers), Depths(layers), &
  !                          COSMIC-parameterset, neutron_integral_AFast, output(cells) )
  !
  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)   :: SoilMoisture" Soil Moisture
  !>        \param[in] "real(dp), dimension(:)   :: Horizons" Horizon depths
  !>        \param[in] "real(dp), dimension(:)   :: neutron_integral_AFast" Tabular for Int Approx
  !>        \param[in] "real(dp)                 :: interc" interception
  !>        \param[in] "real(dp)                 :: snowpack" snowpack
  !>        \param[in] "real(dp)                 :: L1_No_count" L1_No_count
  !>        \param[in] "real(dp), dimension(:)   :: L1_bulkDens" Bulk Density
  !>        \param[in] "real(dp), dimension(:)   :: L1_latticeWater" Lattice Water
  !>        \param[in] "real(dp), dimension(:)   :: L1_COSMICL3" L3 from the COSMIC module

  !
  !     INTENT(INOUT)
  !>        \param[inout] "real(dp)              :: neutrons" Neutron counts
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         Horizons(:) must not be zero.
  !
  !     EXAMPLE
  !         see supplementaries in literature
  !
  !     LITERATURE
  !         J. Shuttleworth, R. Rosolem, M. Zreda, and T. Franz,
  !         The COsmic-ray Soil Moisture Interaction Code (COSMIC) for use in data assimilation,
  !         HESS, 17, 3205-3217, 2013, doi:10.5194/hess-17-3205-2013
  !         Support and Code: http://cosmos.hwr.arizona.edu/Software/cosmic.html
  !
  !     HISTORY
  !>        \author Martin Schroen, originally written by Rafael Rosolem
  !>        \date Mar 2015
  !>         Rohini Kumar                  Oct 2021 - Neutron count module to mHM (5.11.2) - develop branch
  subroutine COSMIC(SoilMoisture, Horizons, neutron_integral_AFast, &
       interc, snowpack, L1_N0, L1_bulkDens, L1_latticeWater, L1_COSMICL3, &
       neutrons )

    use mo_mhm_constants, only: H2Odens, COSMIC_alpha, COSMIC_L1, COSMIC_L2, COSMIC_L4 !, COSMIC_N
    use mo_constants,     only: PI_dp
    implicit none

    real(dp), dimension(:),          intent(in)     :: SoilMoisture
    real(dp), dimension(:),          intent(in)     :: Horizons
    real(dp), dimension(:),          intent(in)     :: neutron_integral_AFast
    real(dp),                        intent(in)     :: interc
    real(dp),                        intent(in)     :: snowpack
    real(dp),                        intent(in)     :: L1_N0
    real(dp), dimension(:),          intent(in)     :: L1_bulkDens
    real(dp), dimension(:),          intent(in)     :: L1_latticeWater
    real(dp), dimension(:),          intent(in)     :: L1_COSMICL3
    real(dp),                        intent(inout)  :: neutrons

    real(dp) :: lambdaHigh
    real(dp) :: lambdaFast
    real(dp) :: totflux
    real(dp) :: sm             ! SoilMoisture
    real(dp) :: lw             ! lattice water
    real(dp) :: bd             ! bulk density
    real(dp) :: L3
    integer(i4):: iFlag_snowlayer_intecept                 ! 1 if snowlayer is active, 0 else

    real(dp), dimension(size(Horizons)+1)   :: zthick      ! Soil layer thickness (cm)
    real(dp), dimension(:), allocatable     :: isoimass    ! Integrated dry soil mass above layer (g)
    real(dp), dimension(:), allocatable     :: iwatmass    ! Integrated water mass above layer (g)
    real(dp), dimension(:), allocatable     :: hiflux      ! High energy neutron flux
    real(dp), dimension(:), allocatable     :: xeff        ! Fast neutron source strength of layer
    real(dp), dimension(:), allocatable     :: h2oeffheight! "Effective" height of water in layer (g/cm3)
    real(dp), dimension(:), allocatable     :: h2oeffdens  ! "Effective" density of water in layer (g/cm3)
    real(dp), dimension(:), allocatable     :: fastflux    ! Contribution to above-ground neutron flux

    integer(i4)                             :: layers      ! Total number of soil layers
    integer(i4)                             :: ll


    ! # of layers
    layers   = size(SoilMoisture) + 1 ! soil horizons + one additional snowpack and interception layer

    ! allocate and initalize
    allocate(hiflux(layers),xeff(layers),&
             h2oeffdens(layers),h2oeffheight(layers),fastflux(layers),&
             isoimass(layers),iwatmass(layers))

    zthick(:)      = 0.0_dp
    isoimass(:)    = 0.0_dp
    iwatmass(:)    = 0.0_dp
    hiflux(:)      = 0.0_dp
    xeff(:)        = 0.0_dp
    h2oeffdens(:)  = 0.0_dp
    h2oeffheight(:)= 0.0_dp
    fastflux(:)    = 0.0_dp
    totflux        = 0.0_dp
    lambdaHigh     = 0.0_dp
    lambdaFast     = 0.0_dp
    sm             = 0.0_dp
    lw             = 0.0_dp
    bd             = 0.0_dp
    L3             = 1.0_dp

    ! switch of snowlayer+Intercept ON and off
    ! ON = 1; OFF = 0
    ! its hardocded here but can be later part of *.nml file
    iFlag_snowlayer_intecept = 0

    ! layer 1 is the surface layer. layer 2 up to layers are the usual layers
    do ll = 1,layers

       ! High energy neutron downward flux
       ! The integration is now performed at the node of each layer (i.e., center of the layer)

       !ToDo: maybe put zthick into global constants, so it is an input paramter
       ! Soil Layers and Thicknesses are constant in mHM, they could be defined outside of this function
       ! except the top layer thickness, which is dependend on the snow for example
       ! zthick will be in cm, as all heigths are in cm in this module
       call layerThickness(ll,Horizons,interc,snowpack,zthick)

       if( ( zthick(ll) .GT. 0.0_dp ) .AND. ( (iFlag_snowlayer_intecept .GT. 0 ) .OR. (ll.NE.1) ) ) then

          call loopConstants(ll, SoilMoisture(:), L1_bulkDens(:), &
               L1_latticeWater(:), L1_COSMICL3(:), sm, bd, lw, L3)

          if (ll .EQ. 1) then
             h2oeffdens(ll) = H2Odens/1000.0_dp
          else
             ! calculate the effective height of water in each layer in cm
             ! because neutron standard measurements are in cm
             call layerWaterHeight(ll,sm,h2oeffheight)
             ! divided by the thickness of the layers,we get the effective density
             h2oeffdens(ll) = (h2oeffheight(ll) + lw/10.0_dp)/zthick(ll)*H2Odens/1000.0_dp
          endif

          ! Assuming an area of 1 cm2
          ! we integrate the bulkdensity/h2oeffdens down to the middle of the layer ll:
          isoimass(ll) = bd*(0.5_dp*zthick(ll))*1.0_dp
          iwatmass(ll) = h2oeffdens(ll)*(0.5_dp*zthick(ll))*1.0_dp
          if ( ll .gt. 1 ) then
            isoimass(ll) = isoimass(ll)+isoimass(ll-1)+bd*(0.5_dp*zthick(ll-1))*1.0_dp
            iwatmass(ll) = iwatmass(ll)+iwatmass(ll-1)+h2oeffdens(ll-1)*(0.5_dp*zthick(ll-1))*1.0_dp
          endif

          lambdaHigh = isoimass(ll)/COSMIC_L1 + iwatmass(ll)/COSMIC_L2
          lambdaFast = isoimass(ll)/L3 + iwatmass(ll)/COSMIC_L4

          hiflux(ll) = exp(-lambdaHigh)
          xeff(ll)   = zthick(ll)*(COSMIC_alpha*bd + h2oeffdens(ll))

          call lookUpIntegral(fastflux(ll),neutron_integral_AFast,lambdaFast)

          ! After contribution from all directions are taken into account,
          ! need to multiply fastflux by 2/pi
          fastflux(ll) = (2.0_dp/PI_dp) * fastflux(ll)

          ! Low energy (fast) neutron upward flux
          totflux = totflux + hiflux(ll) * xeff(ll) * fastflux(ll)
       endif

    enddo

    !  neutrons=COSMIC_N*totflux
    !! >> now based on global parameter given in mhm_paramater.nml
    neutrons = L1_N0*totflux

    !! free space
    deallocate( hiflux, xeff, h2oeffheight, h2oeffdens, fastflux, isoimass, iwatmass)

  end subroutine COSMIC




  ! >>> Loop constants
  subroutine loopConstants(ll, SoilMoisture,L1_bulkDens,L1_latticeWater,&
       L1_COSMICL3,sm,bd,lw,L3 )
     implicit none
     integer(i4),            intent(in)  :: ll
     real(dp), dimension(:), intent(in)  :: SoilMoisture
     real(dp), dimension(:), intent(in)  :: L1_bulkDens
     real(dp), dimension(:), intent(in)  :: L1_latticeWater
     real(dp), dimension(:), intent(in)  :: L1_COSMICL3
     real(dp),               intent(out) :: sm  ! SoilMoisture
     real(dp),               intent(out) :: bd  ! Bulk density
     real(dp),               intent(out) :: lw  ! Lattice water
     real(dp),               intent(out) :: L3  ! L3

     if( ll .EQ. 1 ) then
       !ToDo
       sm = 0.0_dp
       bd = 0.0_dp
       lw = 0.0_dp
       L3 = 1.0_dp
     else
       sm = SoilMoisture(ll-1)
       bd = L1_bulkDens(ll-1)
       lw = L1_latticeWater(ll-1)
       L3 = L1_COSMICL3(ll-1)
     endif
   end subroutine loopConstants



  ! >> layer thickness
  subroutine layerThickness(ll,Horizons,interc,snowpack,zthick)
     implicit none
     integer(i4), intent(in)              :: ll
     real(dp),dimension(:),    intent(in) :: Horizons
     real(dp),                 intent(in) :: interc
     real(dp),                 intent(in) :: snowpack
     real(dp),dimension(:), intent(out)   :: zthick

     if (ll.eq.1) then
        zthick(ll)=(snowpack+interc)/10.0_dp
     else if (ll.eq.2) then
        zthick(ll)=Horizons(ll-1)/10.0_dp
     else
        zthick(ll)=(Horizons(ll-1) - Horizons(ll-2))/10.0_dp
     endif
  end subroutine


  ! >> layer specific water height
  subroutine layerWaterHeight(ll,sm,h2oeffheight)
     implicit none
     integer(i4),           intent(in)  :: ll
     real(dp),              intent(in)  :: sm
     real(dp),dimension(:), intent(out) :: h2oeffheight
    ! The effective water height in each layer in each profile:
    ! ToDo:This should include in future: roots, soil organic matter
    h2oeffheight(ll) = sm/10.0_dp
  end subroutine

  ! integrade a monotonuous function f, dependend on two parameters c and phi
  ! xmin and xmax are the borders for the integration
  ! if the values for f(xmin) or f(xmax) are undefinde (like exp(-1/0)), they
  ! can be set with fxmin, fxmax.
  ! eps is for the accuracy of the result. If the function f is monotonuous, the
  ! error is at most eps.
  ! steps is the maximum number of interpolation points. It is overriding the
  ! error and is the maximum number of steps. A specification of the error
  ! though still has an impact. If the function is interpolated well enough
  ! in a specific flat region regarding the error it can be interpolated better
  ! in a less flat region.
  !
  !For the specific given integral it is very precise with steps=1024
  subroutine approx_mon_int(res,f,c,xmin,xmax,eps,steps,fxmin,fxmax)
     implicit none
     real(dp)                                     :: res
     real(dp), external                           :: f
     real(dp),                        intent(in)  :: c
     real(dp),                        intent(in)  :: xmax
     real(dp),                        intent(in)  :: xmin
     real(dp),                        optional    :: eps
     integer(i4),                     optional    :: steps
     real(dp),                        optional    :: fxmin
     real(dp),                        optional    :: fxmax

     !locale variables
     real(dp)  :: epstemp
     integer(i4) :: stepstemp
     real(dp)  :: fxmintemp
     real(dp)  :: fxmaxtemp

     ! init
     if (.not. present(eps)) then
        epstemp=0.001_dp
     else
        epstemp=eps
     endif

     if (.not. present(steps)) then
        stepstemp=0
     else
        stepstemp=steps
     endif

     if (.not. present(fxmin)) then
        fxmintemp=f(c,xmin)
     else
        fxmintemp=fxmin
     endif

     if (.not. present(fxmax)) then
        fxmaxtemp=f(c,xmax)
     else
        fxmaxtemp=fxmax
     endif

     res=0.0_dp

     if (stepstemp .gt. 0) then
        call approx_mon_int_steps(res,f,c,xmin,xmax,epstemp,stepstemp,fxmintemp,fxmaxtemp)
     else
        call approx_mon_int_eps(res,f,c,xmin,xmax,epstemp,fxmintemp,fxmaxtemp)
     endif

  end subroutine

  recursive subroutine approx_mon_int_steps(res,f,c,xmin,xmax,eps,steps,fxmin,fxmax)
     implicit none
     real(dp)                                     :: res
     real(dp), external                           :: f
     real(dp),                        intent(in)  :: c
     real(dp),                        intent(in)  :: xmax
     real(dp),                        intent(in)  :: xmin
     real(dp),                        intent(in)  :: eps
     integer(i4),                     intent(in)  :: steps
     real(dp),                        intent(in)  :: fxmin
     real(dp),                        intent(in)  :: fxmax

     !locale variables
     real(dp)  :: xm
     real(dp)  :: fxm
     real(dp)  :: err

     xm = (xmax+xmin)/2.0_dp
     fxm= f(c,xm)

     err=abs((fxmax-fxm)*(xmax-xm))
     if ((err .gt. eps).and.(steps .gt. 1)) then
        call approx_mon_int_steps(res,f,c,xm,xmax,eps/2.0,steps-steps/2,fxm,fxmax)
     else
        res=res+(xmax-xm)*(fxmax+fxm)/2.0_dp
     endif

     err=abs((fxm-fxmin)*(xm-xmin))
     if ((err .gt. eps).and.(steps .gt. 1)) then
        call approx_mon_int_steps(res,f,c,xmin,xm,eps/2.0,steps/2,fxmin,fxm)
     else
        res=res+(xm-xmin)*(fxm+fxmin)/2.0_dp
     endif
  end subroutine

  recursive subroutine approx_mon_int_eps(res,f,c,xmin,xmax,eps,fxmin,fxmax)
     implicit none
     real(dp)                                     :: res
     real(dp), external                           :: f
     real(dp),                        intent(in)  :: c
     real(dp),                        intent(in)  :: xmax
     real(dp),                        intent(in)  :: xmin
     real(dp),                        intent(in)  :: eps
     real(dp),                        intent(in)  :: fxmin
     real(dp),                        intent(in)  :: fxmax

     !locale variables
     real(dp)  :: xm
     real(dp)  :: fxm
     real(dp)  :: err

     xm = (xmax+xmin)/2.0_dp
     fxm= f(c,xm)

     err=abs((fxmax-fxm)*(xmax-xm))
     if (err .gt. eps) then
        call approx_mon_int_eps(res,f,c,xm,xmax,eps/2.0,fxm,fxmax)
     else
        res=res+(xmax-xm)*(fxmax+fxm)/2.0_dp
     endif

     err=abs((fxm-fxmin)*(xm-xmin))
     if (err .gt. eps) then
        call approx_mon_int_eps(res,f,c,xmin,xm,eps/2.0,fxmin,fxm)
     else
        res=res+(xm-xmin)*(fxm+fxmin)/2.0_dp
     endif
  end subroutine


  ! -----------------------------------------------------------------------------------
  !     NAME
  !         TabularIntegralAFast
  !
  !     PURPOSE
  !>        \brief Save approximation data for A_fast
  !>        \details The COSMIC subroutine needs A_fast to be calculated.
  !>            A_fast=int_{0}^{pi/2} exp(-Lambda_fast(z)/cos(phi) dphi)
  !>            This subroutine stores data for intsize values for
  !>            c:=Lambda_fast(z) between 0 and maxC, and will be written
  !>            into the global array variable neutron_integral_AFast.
  !>            The calculation of the values is done with a very precise
  !>            recursive approximation subroutine. That recursive subroutine
  !>            should not be used inside the time, cells and layer loops, because
  !>            it is slow.
  !>            Inside the loops in the module COSMIC the tabular is used to
  !>            estimate A_fast, if 0<c<maxC, otherwise the recursive
  !>            approximation is used.
  !         ------------------------------------------------------------------
  !         TabularIntegralAFast: a tabular for calculations with splines
  !         ------------------------------------------------------------------
  !
  !     CALLING SEQUENCE
  !         call TabularIntegralAFast(neutron_integral_AFast,intsize,maxC)
  !
  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:,:) :: SoilMoisture" Soil Moisture
  !>        \param[in] "real(dp), dimension(:)   :: Horizons" Horizon depths
  !>        \param[in] "real(dp), dimension(:)   :: params" ! N0, N1, N2, alpha0, alpha1, L30, L31
  !>        \param[in] "integer(i4)              :: intsize" ! number of values for the approximation
  !>        \param[in] "real(dp)                 :: maxC" ! maximum value for A_fast
  !
  !     INTENT(INOUT)
  !>        \param[out] "real(dp), dimension(intsize) :: neutron_integral_AFast" approximation values
  !
  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         intsize and maxC must be positive
  !
  !     EXAMPLE
  !         intsize=8000, maxC=20.0_dp
  !
  !     LITERATURE
  !         see splines for example
  !
  !     HISTORY
  !>        \author Maren Kaluza
  !>        \date Nov 2017

  subroutine TabularIntegralAFast(integral,maxC)
     use mo_constants, only: PI_dp
     implicit none
     real(dp), dimension(:)              :: integral
     real(dp), intent(in)                :: maxC

     !local variables
     integer(i4)                         :: i
     real(dp)                            :: c
     integer(i4)                         :: intsize

     intsize=size(integral)-2

     do i=1,intsize+1
       c =real(i-1,dp)*maxC/real(intsize,dp)
       call approx_mon_int(integral(i),&
           intgrandFast,c,0.0_dp,PI_dp/2.0_dp,steps=1024,fxmax=0.0_dp)
     enddo
     integral(intsize+2)=maxC
  end subroutine

  ! if c>1.0, the function can be fitted very nice with gnuplot
  ! pi/2*exp(a*x**b)
  subroutine TabularIntegralHermAFast(integral,intsize,maxC)
     use mo_constants, only: PI_dp
     implicit none
     real(dp), dimension(:)              :: integral
     integer(i4), intent(in)             :: intsize
     real(dp), intent(in)                :: maxC

     !local variables
     integer(i4)                         :: i
     real(dp)                            :: c

     do i=1,intsize/2+1
       c =real(i-1,dp)*maxC/real(intsize/2,dp)
       call approx_mon_int(integral(2*i-1),&
           intgrandFast,c,0.0_dp,PI_dp/2.0_dp,steps=1024,fxmax=0.0_dp)
       call approx_mon_int(integral(2*i),&
           intgrandDerivFast,c,0.0_dp,PI_dp/2.0_dp,steps=1024,fxmax=0.0_dp)
           integral(2*i)=integral(2*i)*maxC/real(intsize/2,dp)
     enddo
  end subroutine

  ! if c>1.0, the function can be fitted very nice with gnuplot
  ! pi/2*exp(a*x**b)
  subroutine lookUpIntegral(res,integral,c)
     use mo_constants, only: PI_dp
     implicit none
     real(dp),              intent(out):: res
     real(dp), dimension(:),intent(in) :: integral
     real(dp),              intent(in) :: c

     !local variables
     integer(i4) :: place
     real(dp)    :: mu
     integer(i4) :: intsize
     real(dp)    :: maxC

     intsize=size(integral)-2
     maxC=integral(intsize+2)
     mu=c*real(intsize,dp)/maxC
     place=int(mu,i4)+1
     if (place .gt. intsize) then
       !call approx_mon_int(res,intgrandFast,c,0.0_dp,PI_dp/2.0_dp,steps=1024,fxmax=0.0_dp)
       !write(*,*) 'Warning: Lambda_Fast is huge. Slow integration used.'
       res=(PI_dp/2.0_dp)*exp(-1.57406_dp*c**0.815488_dp)
     else
        mu=mu-real(place-1,dp)
        res=(1.0_dp-mu)*integral(place)+mu*integral(place+1)
        res=res
     end if
  end subroutine

  subroutine lookUpHermiteIntegral(res,integral,intsize,c,maxC)
     use mo_constants, only: PI_dp
     implicit none
     real(dp)                         :: res
     real(dp), dimension(:),intent(in):: integral
     integer(i4), intent(in)          :: intsize
     real(dp), intent(in)             :: c
     real(dp), intent(in)             :: maxC

     !local variables
     integer(i4) :: place
     real(dp)    :: mu

     mu=c*real(intsize/2,dp)/maxC
     place=int(mu,i4)+1
     if (place .gt. intsize) then
       !call approx_mon_int(res,intgrandFast,c,0.0_dp,PI_dp/2.0_dp,steps=1024,fxmax=0.0_dp)
       !write(*,*) 'Warning: Lambda_Fast is huge. Slow integration used.'
       res=(PI_dp/2.0_dp)*exp(-1.57406_dp*c**0.815488_dp)
     else
        mu=mu-real(place-1,dp)
        res=h00(mu)*integral(2*place-1)+h01(mu)*integral(2*place+1)+&
            h10(mu)*integral(2*place  )+h11(mu)*integral(2*place+2)
     end if
  end subroutine

  !very bad approximation for the integral from the COSMIC paper. Maybe there is
  !some copying error? They claim, this function would have an error of less
  !then 1/1000 (which can not be true anyway, because the integral goes to zero
  !for c->\infty, and the last if case is a polynome with some coefficients
  !unequal to zero and therefore tends to \pm \infty), but even in the first 5
  !cases this approximation has sometimes an error of about 1/3 in case 4.
  subroutine COSMICeffIntegration(res,x)
     implicit none
     real(dp)                                     :: res
     real(dp),                        intent(in)  :: x
     !local variables
     real(dp)  :: a,b,c,d

    ! write(*,*) x
     if (x .le. 0.05_dp) then
        a=-347.86105_dp
        b=41.64233_dp
        c=-4.018_dp
        d=-0.00018_dp
        res=expPolynomDeg3(x,a,b,c,d)
      !  write(*,*) 'c1'
     else if (x .le. 0.1_dp) then
        a=-16.24066_dp
        b=6.64468_dp
        c=-2.82003_dp
        d=-0.01389_dp
        res=expPolynomDeg3(x,a,b,c,d)
      !  write(*,*) 'c2'
     else if (x .le. 0.5_dp) then
        a=-0.95245_dp
        b=1.44751_dp
        c=-2.18933_dp
        d=-0.04034_dp
        res=expPolynomDeg3(x,a,b,c,d)
      !  write(*,*) 'c3'
     else if (x .le. 1.0_dp) then
        a=-0.09781_dp
        b=0.36907_dp
        c=-1.72912_dp
        d=-0.10761_dp
        res=expPolynomDeg3(x,a,b,c,d)
      !  write(*,*) 'c4'
     else if (x .le. 5.0_dp) then
        a=-0.00416_dp
        b=0.05808_dp
        c=-1.361482_dp
        d=-0.25822_dp
        res=expPolynomDeg3(x,a,b,c,d)
      !  write(*,*) 'c5'
     else
        a=0.0_dp
        b=0.00061_dp
        c=-1.04847_dp
        d=-0.96617_dp
        res=expPolynomDeg3(x,a,b,c,d)
      !  write(*,*) 'c6'
     endif

  end subroutine

  subroutine oldIntegration(res,c)
     use mo_constants, only: PI_dp
     implicit none
     real(dp)                                     :: res
     real(dp),                        intent(in)  :: c

     ! local variables
     real(dp) :: zdeg
     real(dp) :: zrad
     real(dp) :: costheta
     real(dp) :: dtheta

     integer(i4) :: angle ! loop indices for an integration interval
     ! Angle distribution parameters (HARDWIRED)
     ! rr: Using 0.5 deg angle intervals appears to be sufficient
     ! rr: (smaller angles increase the computing time for COSMIC)

     dtheta = 0.5_dp*(PI_dp/180.0_dp)

     ! This second loop needs to be done for the distribution of angles for fast neutron release
     ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
     ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.
     res = 0.0_dp
     do angle=0,179
        zdeg     = real(angle,dp)*0.5_dp
        zrad     = (zdeg*PI_dp)/180.0_dp
        costheta = cos(zrad)
        ! Angle-dependent low energy (fast) neutron upward flux
        res  = res + exp(-c/costheta)*dtheta
     enddo
  end subroutine

  function intgrandFast(c,phi)
     implicit none
     real(dp) :: intgrandFast
     real(dp), intent(in) :: c
     real(dp), intent(in) :: phi
     intgrandFast=exp(-c/cos(phi))
     return
  end function

  function intgrandDerivFast(c,phi)
     implicit none
     real(dp) :: intgrandDerivFast
     real(dp), intent(in) :: c
     real(dp), intent(in) :: phi
     intgrandDerivFast=(-1.0_dp/cos(phi))*exp(-c/cos(phi))
     return
  end function

  function expPolynomDeg3(x,a,b,c,d)
     implicit none
     real(dp) :: expPolynomDeg3
     real(dp), intent(in) :: x
     real(dp), intent(in) :: a,b,c,d

     expPolynomDeg3=exp(a*x**3+b*x**2+c*x+d)
     return
  end function

  ! hermite polynoms
  function h00(t)
     implicit none
     real(dp)             :: h00
     real(dp), intent(in) :: t
     h00=2.0_dp*t**3-3.0_dp*t**2+1.0_dp
     return
  end function

  function h01(t)
     implicit none
     real(dp)             :: h01
     real(dp), intent(in) :: t
     h01=-2.0_dp*t**3+3.0_dp*t**2
     return
  end function

  function h10(t)
     implicit none
     real(dp)             :: h10
     real(dp), intent(in) :: t
     h10=t**3-2.0_dp*t**2+t
     return
  end function

  function h11(t)
     implicit none
     real(dp)             :: h11
     real(dp), intent(in) :: t
     h11=t**3-t**2
     return
  end function

END MODULE mo_neutrons
