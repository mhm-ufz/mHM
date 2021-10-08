!>       \file mo_neutrons.f90

!>       \brief Models to predict neutron intensities above soils

!>       \details The number of neutrons above the ground is directly related to
!>       the number soil water content in the ground, air, vegetation and/or snow.
!>       This module forward-models neutron abundance as a state variable for each cell.

!>       \authors Martin Schroen

!>       \date Mar 2015
!>       THIS MODULE IS WORK IN PROGRESS, DO NOT USE FOR RESEARCH.

! Modifications:

MODULE mo_neutrons

  ! Written  Martin Schroen, Mar 2015
  ! Modified 

  USE mo_kind, ONLY : i4, dp
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
  !    NAME
  !        DesiletsN0

  !    PURPOSE
  !>       \brief Calculate neutrons from soil moisture in the first layer.

  !>       \details Using the N0-relation derived by Desilets, neutron
  !>       counts above the ground (one value per cell in mHM) can be
  !>       derived by a semi-empirical, semi-physical relation.
  !>       The result depends on N0, the neutron counts for 0% soil mositure.
  !>       This variable is site-specific and is a global parameter in mHM.

  !>       N0 formula based on Desilets et al. 2010

  !>       Horizons(1) must not be zero.

  !>       N0=1500cph, SoilMoisture(1,1)=700mm, Horizons(1)=200mm
  !>       1500*(0.372+0.0808/ (70mm/200mm + 0.115))
  !>       DesiletsN0 = 819cph

  !>       Desilets, D., M. Zreda, and T. P. A. Ferre (2010),
  !>       Nature's neutron probe: Land surface hydrology at an elusive scale
  !>       with cosmic rays, WRR, 46, W11505, doi:10.1029/2009WR008726.

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: SoilMoisture" Soil Moisture
  !>       \param[in] "real(dp), dimension(:) :: Horizons"     Horizon depths
  !>       \param[in] "real(dp) :: N0"                         dry neutron counts

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp) :: neutrons" Neutron counts

  !    HISTORY
  !>       \authors Martin Schroen

  !>       \date Mar 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine DesiletsN0(SoilMoisture, Horizons, N0, neutrons)

    use mo_mhm_constants, only : Desilets_a0, Desilets_a1, Desilets_a2

    implicit none

    ! Soil Moisture
    real(dp), dimension(:), intent(in) :: SoilMoisture

    ! Horizon depths
    real(dp), dimension(:), intent(in) :: Horizons

    ! dry neutron counts
    real(dp), intent(in) :: N0

    ! Neutron counts
    real(dp), intent(inout) :: neutrons


    ! only use first soil layer
    neutrons = N0 * (Desilets_a1 + Desilets_a0 / (SoilMoisture(1) / Horizons(1) + Desilets_a2))

  end subroutine DesiletsN0

  ! -----------------------------------------------------------------------------------
  !    NAME
  !        COSMIC

  !    PURPOSE
  !>       \brief Calculate neutrons from soil moisture in all layers.

  !>       \details Neutron counts above the ground (one value per cell in mHM)
  !>       can be derived by a simplified physical neutron transport simulation.
  !>       Fast cosmic-Ray neutrons are generated in the soil and attenuated
  !>       differently in water and soil. The remaining neutrons that reached
  !>       the surface relate to the profile of soil water content below.
  !>       Variables like N, alpha and L3 are site-specific and need to be calibrated.
  !>       ADDITIONAL INFORMATION

  !>       COSMIC model based on Shuttleworth et al. 2013

  !>       Horizons(:) must not be zero.

  !>       see supplementaries in literature

  !>       J. Shuttleworth, R. Rosolem, M. Zreda, and T. Franz,
  !>       The COsmic-ray Soil Moisture Interaction Code (COSMIC) for use in data assimilation,
  !>       HESS, 17, 3205-3217, 2013, doi:10.5194/hess-17-3205-2013
  !>       Support and Code: http://cosmos.hwr.arizona.edu/Software/cosmic.html

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: SoilMoisture"           Soil Moisture
  !>       \param[in] "real(dp), dimension(:) :: Horizons"               Horizon depths
  !>       \param[in] "real(dp), dimension(:) :: params"                 ! N0, N1, N2, alpha0, alpha1, L30, L31
  !>       \param[in] "real(dp), dimension(:) :: neutron_integral_AFast" Tabular for Int Approx

  !    INTENT(INOUT)
  !>       \param[inout] "real(dp) :: neutrons" Neutron counts

  !    HISTORY
  !>       \authors Martin Schroen, originally written by Rafael Rosolem

  !>       \date Mar 2015

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine COSMIC(SoilMoisture, Horizons, params, neutron_integral_AFast, &
                       L1_bulkDens, &
                       L1_latticeWater, &
                       L1_COSMICL3, &
                       interc              , & ! Interception
                       snowpack            , & ! Snowpack
                       neutrons)
    
    use mo_constants, only: PI_dp
	use mo_mhm_constants, only: H2Odens, &
        COSMIC_N, COSMIC_alpha, COSMIC_bd, COSMIC_vwclat, &
        COSMIC_L1, COSMIC_L2, COSMIC_L3, COSMIC_L4
		

    implicit none

    ! Soil Moisture
    real(dp), dimension(:), intent(in) :: SoilMoisture

    ! Horizon depths
    real(dp), dimension(:), intent(in) :: Horizons

    ! ! N0, N1, N2, alpha0, alpha1, L30, L31
    real(dp), dimension(:), intent(in) :: params

    ! Tabular for Int Approx
    real(dp), dimension(:), intent(in) :: neutron_integral_AFast
	
	real(dp), dimension(:),          intent(in)     :: L1_bulkDens
    real(dp), dimension(:),          intent(in)     :: L1_latticeWater
    real(dp), dimension(:),          intent(in)     :: L1_COSMICL3
    real(dp),                        intent(in)     :: interc
    real(dp),                        intent(in)     :: snowpack

    ! Neutron counts
    real(dp), intent(inout) :: neutrons

    real(dp) :: lambdaHigh

    real(dp) :: lambdaFast

    real(dp) :: totflux
	real(dp) :: sm             ! SoilMoisture
    real(dp) :: lw             ! lattice water
    real(dp) :: bd             ! bulk density
    real(dp) :: L3
    integer(i4):: snowlayer    ! 1 if snowlayer is active, 0 else


    ! Soil layer thickness (cm)
    real(dp), dimension(size(Horizons)) :: zthick

    ! Integrated dry soil mass above layer (g)
    real(dp), dimension(:), allocatable :: isoimass

    ! Integrated water mass above layer (g)
    real(dp), dimension(:), allocatable :: iwatmass

    ! High energy neutron flux
    real(dp), dimension(:), allocatable :: hiflux

    ! Fast neutron source strength of layer
    real(dp), dimension(:), allocatable :: xeff
	
	! "Effective" height of water in layer (g/cm3)
	real(dp), dimension(:), allocatable     :: h2oeffheight

    ! "Effective" density of water in layer (g/cm3)
    real(dp), dimension(:), allocatable :: h2oeffdens

    ! Contribution to above-ground neutron flux
    real(dp), dimension(:), allocatable :: fastflux

    ! Total number of soil layers
    integer(i4) :: layers = 1

    integer(i4) :: ll = 1


   layers   = size(SoilMoisture)+1 ! 3, one additional snowpack layer
    
    allocate(hiflux(layers),xeff(layers),&
             h2oeffdens(layers),h2oeffheight(layers),fastflux(layers),&
             isoimass(layers),iwatmass(layers))

    zthick(:) = 0.0_dp * params(1) ! <-- this multiplication with params(1) is not needed, only to make params USED
    !                                     !     PLEASE remove when possible 
    isoimass(:) = 0.0_dp
    iwatmass(:) = 0.0_dp
    hiflux(:) = 0.0_dp
    xeff(:) = 0.0_dp
    h2oeffdens(:) = 0.0_dp
	h2oeffheight(:)= 0.0_dp
    fastflux(:) = 0.0_dp
    lambdaHigh = 0.0_dp
    lambdaFast = 0.0_dp
    totflux = 0.0_dp
	sm = 0.0_dp
    lw = 0.0_dp
    bd = 0.0_dp
    L3 = 1.0_dp

     snowlayer=0
    ! Soil Layers and Thicknesses are constant in mHM, they could be defined outside of this function
    zthick(1) = Horizons(1) / 10.0_dp - 0.0_dp
    do ll = 2, layers
      zthick(ll) = (Horizons(ll) - Horizons(ll - 1)) / 10.0_dp
    enddo

    do ll = 1, layers

       ! High energy neutron downward flux
       ! The integration is now performed at the node of each layer (i.e., center of the layer)

       !ToDo: maybe put zthick into global constants, so it is an input paramter
       ! Soil Layers and Thicknesses are constant in mHM, they could be defined outside of this function
       ! except the top layer thickness, which is dependend on the snow for example
       ! zthick will be in cm, as all heigths are in cm in this module
       call layerThickness(ll,Horizons,interc,snowpack,zthick)

       if (zthick(ll).gt.0.0_dp .and. (snowlayer.gt.0 .or. ll.ne.1)) then
          call loopConstants(ll,&
                    SoilMoisture(:),L1_bulkDens(:),L1_latticeWater(:),&
                    L1_COSMICL3(:),sm,bd,lw,L3)


          if (ll.eq.1) then
             h2oeffdens(ll) = H2Odens/1000.0_dp
          else
             ! calculate the effective height of water in each layer in cm
             ! because neutron standard measurements are in cm
             call layerWaterHeight(ll,sm,h2oeffheight)
             ! divided by the thickness of the layers,we get the effective density
             h2oeffdens(ll) = (h2oeffheight(ll) +lw/10.0_dp)*H2Odens/zthick(ll)/1000.0_dp  
          endif

          ! Assuming an area of 1 cm2
          ! we integrate the bulkdensity/h2oeffdens down to the middle of the layer ll:
          isoimass(ll) = bd*(0.5_dp*zthick(ll))*1.0_dp 
          iwatmass(ll) = h2oeffdens(ll)*(0.5_dp*zthick(ll))*1.0_dp
          if (ll.gt.1) then
            isoimass(ll) = isoimass(ll)+isoimass(ll-1)+bd*(0.5_dp*zthick(ll-1))*1.0_dp
            iwatmass(ll) = iwatmass(ll)+iwatmass(ll-1)+h2oeffdens(ll-1)*(0.5_dp*zthick(ll-1))*1.0_dp
          endif

          lambdaHigh = isoimass(ll)/COSMIC_L1 + iwatmass(ll)/COSMIC_L2
          lambdaFast = isoimass(ll)/L3 + iwatmass(ll)/COSMIC_L4

          hiflux(ll)  = exp(-lambdaHigh)
          xeff(ll) = zthick(ll)*(COSMIC_alpha*bd + h2oeffdens(ll))

          call lookUpIntegral(fastflux(ll),neutron_integral_AFast,lambdaFast)

          ! After contribution from all directions are taken into account,
          ! need to multiply fastflux by 2/pi
          fastflux(ll)=(2.0_dp/PI_dp)*fastflux(ll)

          ! Low energy (fast) neutron upward flux
          totflux=totflux+hiflux(ll)*xeff(ll)*fastflux(ll)

       endif
    enddo
	
    neutrons = COSMIC_N * totflux

    deallocate( hiflux,&
           xeff, h2oeffheight, h2oeffdens, fastflux,&
           isoimass, iwatmass)

  end subroutine COSMIC
  
  subroutine loopConstants(ll,&
                    SoilMoisture,L1_bulkDens,L1_latticeWater,&
                    L1_COSMICL3,sm,bd,lw,L3)
     implicit none
     integer(i4), intent(in)                    :: ll
     real(dp), dimension(:),        intent(in)  :: SoilMoisture
     real(dp), dimension(:),        intent(in)  :: L1_bulkDens
     real(dp), dimension(:),        intent(in)  :: L1_latticeWater
     real(dp), dimension(:),        intent(in)  :: L1_COSMICL3
     real(dp) :: sm  ! SoilMoisture
     real(dp) :: bd  ! bulk density
     real(dp) :: lw  ! lattice water
     real(dp) :: L3

     if (ll.eq.1) then
       !ToDo
       sm=0.0_dp
       bd=0.0_dp
       lw=0.0_dp
       L3=1.0_dp
     else
       sm=SoilMoisture(ll-1)
       bd=L1_bulkDens(ll-1)
       lw=L1_latticeWater(ll-1)
       L3=L1_COSMICL3(ll-1)
     endif
  end subroutine

  subroutine layerThickness(ll,Horizons,interc,snowpack,zthick)
     implicit none
     integer(i4), intent(in)              :: ll
     real(dp),dimension(:),    intent(in) :: Horizons
     real(dp),                 intent(in) :: interc
     real(dp),                 intent(in) :: snowpack
     real(dp),dimension(:)                :: zthick
       if (ll.eq.1) then
          zthick(ll)=(snowpack+interc)/10.0_dp
       else if (ll.eq.2) then
          zthick(ll)=Horizons(ll-1)/10.0_dp
       else
          zthick(ll)=(Horizons(ll-1)-Horizons(ll-2))/10.0_dp
       endif
  end subroutine

  subroutine layerWaterHeight(ll,sm,h2oeffheight)
     implicit none
     integer(i4), intent(in) :: ll
     real(dp),    intent(in) :: sm
     real(dp),dimension(:)   :: h2oeffheight
    ! The effective water height in each layer in each profile:
    ! ToDo:This should include in future: roots, soil organic matter 
    h2oeffheight(ll) = sm/10.0_dp
  end subroutine


  !    NAME
  !        oldIntegration

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: c"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine oldIntegration(res, c)

    use mo_constants, only : PI_dp

    implicit none

    real(dp) :: res

    real(dp), intent(in) :: c

    real(dp) :: zdeg

    real(dp) :: zrad

    real(dp) :: costheta

    real(dp) :: dtheta

    ! loop indices for an integration interval
    integer(i4) :: angle


    ! Angle distribution parameters (HARDWIRED)
    ! rr: Using 0.5 deg angle intervals appears to be sufficient
    ! rr: (smaller angles increase the computing time for COSMIC)
    dtheta = 0.5_dp * (PI_dp / 180.0_dp)

    ! This second loop needs to be done for the distribution of angles for fast neutron release
    ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
    ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.

    res = 0.0_dp
    do angle = 0, 179
      zdeg = real(angle, dp) * 0.5_dp
      zrad = (zdeg * PI_dp) / 180.0_dp
      costheta = cos(zrad)

      ! Angle-dependent low energy (fast) neutron upward flux
      res = res + exp(-c / costheta) * dtheta
    enddo
  end subroutine

  ! -----------------------------------------------------------------------------------
  !    NAME
  !        TabularIntegralAFast

  !    PURPOSE
  !>       \brief Save approximation data for A_fast

  !>       \details The COSMIC subroutine needs A_fast to be calculated.
  !>       A_fast=int_{0}^{pi/2} exp(-Lambda_fast(z)/cos(phi) dphi)
  !>       This subroutine stores data for intsize values for
  !>       c=Lambda_fast(z) between 0 and maxC, and will be written
  !>       into the global array variable neutron_integral_AFast.
  !>       The calculation of the values is done with a very precise
  !>       recursive approximation subroutine. That recursive subroutine
  !>       should not be used inside the time, cells and layer loops, because
  !>       it is slow.
  !>       Inside the loops in the module COSMIC the tabular is used to
  !>       estimate A_fast, if 0<c<maxC, otherwise the recursive
  !>       approximation is used.

  !>       TabularIntegralAFast: a tabular for calculations with splines

  !>       intsize and maxC must be positive

  !>       intsize=8000, maxC=20.0_dp

  !>       see splines for example

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: maxC" ! maximum value for A_fast

  !    HISTORY
  !>       \authors Maren Kaluza

  !>       \date Nov 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine TabularIntegralAFast(integral, maxC)

    use mo_constants, only : PI_dp

    implicit none

    real(dp), dimension(:) :: integral

    ! ! maximum value for A_fast
    real(dp), intent(in) :: maxC

    integer(i4) :: i

    real(dp) :: c

    integer(i4) :: intsize


    intsize = size(integral) - 2

    do i = 1, intsize + 1
      c = real(i - 1, dp) * maxC / real(intsize, dp)
      call approx_mon_int(integral(i), &
              intgrandFast, c, 0.0_dp, PI_dp / 2.0_dp, steps = 1024, fxmax = 0.0_dp)
    enddo
    integral(intsize + 2) = maxC
  end subroutine

  !For the specific given integral it is very precise with steps=1024
  !    NAME
  !        approx_mon_int

  !    PURPOSE
  !>       \brief TODO: add description
  ! integrate a monotonuous function f, dependend on two parameters c and phi
  ! xmin and xmax are the borders for the integration
  ! if the values for f(xmin) or f(xmax) are undefined (like exp(-1/0)), they
  ! can be set with fxmin, fxmax.
  ! eps is for the accuracy of the result. If the function f is monotonuous, the
  ! error is at most eps.
  ! steps is the maximum number of interpolation points. It is overriding the
  ! error and is the maximum number of steps. A specification of the error
  ! though still has an impact. If the function is interpolated well enough
  ! in a specific flat region regarding the error it can be interpolated better
  ! in a less flat region.
  !

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: c"
  !>       \param[in] "real(dp) :: xmin"
  !>       \param[in] "real(dp) :: xmax"

  !    INTENT(IN), OPTIONAL
  !>       \param[in] "real(dp), optional :: eps"
  !>       \param[in] "integer(i4), optional :: steps"
  !>       \param[in] "real(dp), optional :: fxmin"
  !>       \param[in] "real(dp), optional :: fxmax"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine approx_mon_int(res, f, c, xmin, xmax, eps, steps, fxmin, fxmax)
    implicit none

    real(dp) :: res

    real(dp), external :: f

    real(dp), intent(in) :: c

    real(dp), intent(in) :: xmax

    real(dp), intent(in) :: xmin

    real(dp), intent(in), optional :: eps

    integer(i4), intent(in), optional :: steps

    real(dp), intent(in), optional :: fxmin

    real(dp), intent(in), optional :: fxmax

    real(dp) :: epstemp

    integer(i4) :: stepstemp

    real(dp) :: fxmintemp

    real(dp) :: fxmaxtemp


    ! init
    if (.not. present(eps)) then
      epstemp = 0.001_dp
    else
      epstemp = eps
    endif

    if (.not. present(steps)) then
      stepstemp = 0
    else
      stepstemp = steps
    endif

    if (.not. present(fxmin)) then
      fxmintemp = f(c, xmin)
    else
      fxmintemp = fxmin
    endif

    if (.not. present(fxmax)) then
      fxmaxtemp = f(c, xmax)
    else
      fxmaxtemp = fxmax
    endif

    res = 0.0_dp

    if (stepstemp .gt. 0) then
      call approx_mon_int_steps(res, f, c, xmin, xmax, epstemp, stepstemp, fxmintemp, fxmaxtemp)
    else
      call approx_mon_int_eps(res, f, c, xmin, xmax, epstemp, fxmintemp, fxmaxtemp)
    endif

  end subroutine

  !    NAME
  !        approx_mon_int_steps

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: c"
  !>       \param[in] "real(dp) :: xmin"
  !>       \param[in] "real(dp) :: xmax"
  !>       \param[in] "real(dp) :: eps"
  !>       \param[in] "integer(i4) :: steps"
  !>       \param[in] "real(dp) :: fxmin"
  !>       \param[in] "real(dp) :: fxmax"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  recursive subroutine approx_mon_int_steps(res, f, c, xmin, xmax, eps, steps, fxmin, fxmax)
    implicit none

    real(dp) :: res

    real(dp), external :: f

    real(dp), intent(in) :: c

    real(dp), intent(in) :: xmax

    real(dp), intent(in) :: xmin

    real(dp), intent(in) :: eps

    integer(i4), intent(in) :: steps

    real(dp), intent(in) :: fxmin

    real(dp), intent(in) :: fxmax

    real(dp) :: xm

    real(dp) :: fxm

    real(dp) :: err


    xm = (xmax + xmin) / 2.0_dp
    fxm = f(c, xm)

    err = abs((fxmax - fxm) * (xmax - xm))
    if ((err .gt. eps).and.(steps .gt. 1)) then
      call approx_mon_int_steps(res, f, c, xm, xmax, eps / 2.0, steps - steps / 2, fxm, fxmax)
    else
      res = res + (xmax - xm) * (fxmax + fxm) / 2.0_dp
    endif

    err = abs((fxm - fxmin) * (xm - xmin))
    if ((err .gt. eps).and.(steps .gt. 1)) then
      call approx_mon_int_steps(res, f, c, xmin, xm, eps / 2.0, steps / 2, fxmin, fxm)
    else
      res = res + (xm - xmin) * (fxm + fxmin) / 2.0_dp
    endif
  end subroutine

  !    NAME
  !        approx_mon_int_eps

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: c"
  !>       \param[in] "real(dp) :: xmin"
  !>       \param[in] "real(dp) :: xmax"
  !>       \param[in] "real(dp) :: eps"
  !>       \param[in] "real(dp) :: fxmin"
  !>       \param[in] "real(dp) :: fxmax"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  recursive subroutine approx_mon_int_eps(res, f, c, xmin, xmax, eps, fxmin, fxmax)
    implicit none

    real(dp) :: res

    real(dp), external :: f

    real(dp), intent(in) :: c

    real(dp), intent(in) :: xmax

    real(dp), intent(in) :: xmin

    real(dp), intent(in) :: eps

    real(dp), intent(in) :: fxmin

    real(dp), intent(in) :: fxmax

    real(dp) :: xm

    real(dp) :: fxm

    real(dp) :: err


    xm = (xmax + xmin) / 2.0_dp
    fxm = f(c, xm)

    err = abs((fxmax - fxm) * (xmax - xm))
    if (err .gt. eps) then
      call approx_mon_int_eps(res, f, c, xm, xmax, eps / 2.0, fxm, fxmax)
    else
      res = res + (xmax - xm) * (fxmax + fxm) / 2.0_dp
    endif

    err = abs((fxm - fxmin) * (xm - xmin))
    if (err .gt. eps) then
      call approx_mon_int_eps(res, f, c, xmin, xm, eps / 2.0, fxmin, fxm)
    else
      res = res + (xm - xmin) * (fxm + fxmin) / 2.0_dp
    endif
  end subroutine

  !    NAME
  !        lookUpIntegral

  !    PURPOSE
  !>       \brief TODO: add description
  ! if c>1.0, the function can be fitted very nice with gnuplot
  ! pi/2*exp(a*x**b)

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(:) :: integral"
  !>       \param[in] "real(dp) :: c"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine lookUpIntegral(res, integral, c)

    use mo_constants, only : PI_dp

    implicit none

    real(dp) :: res

    real(dp), dimension(:), intent(in) :: integral

    real(dp), intent(in) :: c

    integer(i4) :: place

    real(dp) :: mu

    integer(i4) :: intsize

    real(dp) :: maxC


    intsize = size(integral) - 2
    maxC = integral(intsize + 2)
    mu = c * real(intsize, dp) / maxC
    place = int(mu, i4) + 1
    if (place .gt. intsize) then
      !call approx_mon_int(res,intgrandFast,c,0.0_dp,PI_dp/2.0_dp,steps=1024,fxmax=0.0_dp)
      !write(*,*) 'Warning: Lambda_Fast is huge. Slow integration used.'
      res = (PI_dp / 2.0_dp) * exp(-1.57406_dp * c**0.815488_dp)
    else
      mu = mu - real(place - 1, dp)
      res = (1.0_dp - mu) * integral(place) + mu * integral(place + 1)
      res = res
    end if
  end subroutine

  !    NAME
  !        intgrandFast

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "real(dp) :: c"
  !>       \param[in] "real(dp) :: phi"

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Jun 2018

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  function intgrandFast(c, phi)
    implicit none

    real(dp) :: intgrandFast

    real(dp), intent(in) :: c

    real(dp), intent(in) :: phi


    intgrandFast = exp(-c / cos(phi))
    return
  end function

END MODULE mo_neutrons
