!> \file mo_neutrons.f90
!> \brief Models to predict neutron intensities above soils

!> \details The number of neutrons above the ground is directly related to
!> the number soil water content in the ground, air, vegetation and/or snow.
!> This module forward-models neutron abundance as a state variable for each cell.

!> \authors Martin Schroen
!> \date Mar 2015

!> THIS MODULE IS WORK IN PROGRESS, DO NOT USE FOR RESEARCH.

! TODO make it faster with pre-calculated horizons and variables
! TODO use global parameters as linear model 

MODULE mo_neutrons

  ! Written  Martin Schroen, Mar 2015
  ! Modified 
  
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
  !>        \brief Calculate neutrons from soil moisture in the first layer.
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
  !>        \param[in] "real(dp), dimension(:,:) :: SoilMoisture" Soil Moisture
  !>        \param[in] "real(dp), dimension(:)   :: Horizons" Horizon depths
  !>        \param[in] "real(dp)                 :: N0" dry neutron counts
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
  !         N0=1500cph, SoilMoisture(1,1)=700mm, Horizons(1)=200mm
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

  subroutine DesiletsN0(SoilMoisture, Horizons, N0, neutrons)

    use mo_mhm_constants, only: Desilets_a0, Desilets_a1, Desilets_a2
    implicit none
    
    real(dp), dimension(:),          intent(in)  :: SoilMoisture
    real(dp), dimension(:),          intent(in)  :: Horizons
    real(dp),                        intent(in)  :: N0          ! from global parameters
    real(dp),                        intent(inout) :: neutrons
    
    ! only use first soil layer
    neutrons = N0 * ( Desilets_a1 + Desilets_a0 / (SoilMoisture(1)/Horizons(1) + Desilets_a2))
  
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
  !>        \param[in] "real(dp), dimension(:,:) :: SoilMoisture" Soil Moisture
  !>        \param[in] "real(dp), dimension(:)   :: Horizons" Horizon depths
  !>        \param[in] "real(dp), dimension(:)   :: params" ! N0, N1, N2, alpha0, alpha1, L30, L31
  !>        \param[in] "real(dp), dimension(:)   :: neutron_integral_AFast" Tabular for Int Approx
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
  
  subroutine COSMIC(SoilMoisture, Horizons, params, neutron_integral_AFast, neutrons)
    
    use mo_mhm_constants, only: H2Odens, &
        COSMIC_bd, COSMIC_vwclat, COSMIC_N, COSMIC_alpha, &
        COSMIC_L1, COSMIC_L2, COSMIC_L3, COSMIC_L4
    use mo_constants, only: PI_dp
    implicit none
    
    real(dp), dimension(:),          intent(in)  :: SoilMoisture
    real(dp), dimension(:),          intent(in)  :: Horizons
    real(dp), dimension(:),          intent(in)  :: params ! 1: N0, 2: N1, 3: N2, 4: alpha0, 5: alpha1, 6: L30, 7. L31
    real(dp), dimension(:),          intent(in)  :: neutron_integral_AFast
    real(dp),                        intent(inout) :: neutrons

    ! local variables
    real(dp) :: lambdaHigh
    real(dp) :: lambdaFast
    real(dp) :: totflux
   
    real(dp), dimension(size(Horizons))   :: zthick      ! Soil layer thickness (cm)
    real(dp), dimension(:), allocatable   :: isoimass    ! Integrated dry soil mass above layer (g)
    real(dp), dimension(:), allocatable   :: iwatmass    ! Integrated water mass above layer (g)
    real(dp), dimension(:), allocatable   :: hiflux      ! High energy neutron flux
    real(dp), dimension(:), allocatable   :: xeff        ! Fast neutron source strength of layer
    real(dp), dimension(:), allocatable   :: h2oeffdens  ! "Effective" density of water in layer (g/cm3)
    real(dp), dimension(:), allocatable   :: fastflux    ! Contribution to above-ground neutron flux

    !
    integer(i4) :: layers=1                 ! Total number of soil layers
    integer(i4) :: ll=1
    !
    layers   = size(SoilMoisture) ! 2
    
    allocate(hiflux(layers),xeff(layers),&
             h2oeffdens(layers),&
             fastflux(layers),&
             isoimass(layers),iwatmass(layers))

    zthick(:)        = 0.0_dp * params(1) ! <-- this multiplication with params(1) is not needed, only to make params USED
    !                                     !     PLEASE remove when possible 
    isoimass(:)      = 0.0_dp
    iwatmass(:)      = 0.0_dp
    hiflux(:)        = 0.0_dp
    xeff(:)          = 0.0_dp
    h2oeffdens(:)    = 0.0_dp
    fastflux(:)      = 0.0_dp
    lambdaHigh       = 0.0_dp
    lambdaFast       = 0.0_dp
    totflux          = 0.0_dp
    
    ! Soil Layers and Thicknesses are constant in mHM, they could be defined outside of this function
    zthick(1) = Horizons(1)/10.0_dp - 0.0_dp
    do ll = 2,layers
       zthick(ll) = (Horizons(ll) - Horizons(ll-1))/10.0_dp
    enddo

    do ll = 1,layers
          
       ! High energy neutron downward flux
       ! The integration is now performed at the node of each layer (i.e., center of the layer)
       h2oeffdens(ll) = ((SoilMoisture(ll) / zthick(ll) / 10.0_dp +COSMIC_vwclat)*H2Odens)/1000.0_dp  

       ! Assuming an area of 1 cm2
       if(ll > 1) then
          isoimass(ll) = isoimass(ll-1) + COSMIC_bd*(0.5_dp*zthick(ll-1))*1.0_dp &
                                          + COSMIC_bd*(0.5_dp*zthick(ll))*1.0_dp
          iwatmass(ll) = iwatmass(ll-1) + h2oeffdens(ll-1)*(0.5_dp*zthick(ll-1))*1.0_dp &
                                          + h2oeffdens(ll)*(0.5_dp*zthick(ll))*1.0_dp
       else
          isoimass(ll) = COSMIC_bd*(0.5_dp*zthick(ll))*1.0_dp 
          iwatmass(ll) = h2oeffdens(ll)*(0.5_dp*zthick(ll))*1.0_dp
       end if

       lambdaHigh = isoimass(ll)/COSMIC_L1 + iwatmass(ll)/COSMIC_L2
       lambdaFast = isoimass(ll)/COSMIC_L3 + iwatmass(ll)/COSMIC_L4
       
       hiflux(ll)  = exp(-lambdaHigh)
       xeff(ll) = zthick(ll)*(COSMIC_alpha*COSMIC_bd + h2oeffdens(ll))
       
       call lookUpIntegral(fastflux(ll),neutron_integral_AFast,lambdaFast)

       ! After contribution from all directions are taken into account,
       ! need to multiply fastflux by 2/pi
       fastflux(ll) = (2.0_dp/PI_dp)*fastflux(ll)

       ! Low energy (fast) neutron upward flux
       totflux = totflux + hiflux(ll)*xeff(ll)*fastflux(ll)

    enddo
    neutrons=COSMIC_N*totflux
    

    deallocate(hiflux,&
           xeff, h2oeffdens, fastflux,&
           isoimass, iwatmass)
           
  end subroutine COSMIC

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

     dtheta   = 0.5_dp*(PI_dp/180.0_dp)

     ! This second loop needs to be done for the distribution of angles for fast neutron release
     ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
     ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.  

     res=0.0_dp
     do angle=0,179
        zdeg     = real(angle,dp)*0.5_dp
        zrad     = (zdeg*PI_dp)/180.0_dp
        costheta = cos(zrad)

        ! Angle-dependent low energy (fast) neutron upward flux
        res  = res + exp(-c/costheta)*dtheta
     enddo
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
  !>            c=Lambda_fast(z) between 0 and maxC, and will be written
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

  ! if c>1.0, the function can be fitted very nice with gnuplot
  ! pi/2*exp(a*x**b)
  subroutine lookUpIntegral(res,integral,c)
     use mo_constants, only: PI_dp
     implicit none
     real(dp)                         :: res
     real(dp), dimension(:),intent(in):: integral
     real(dp), intent(in)             :: c

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

  function intgrandFast(c,phi)
     implicit none
     real(dp) :: intgrandFast
     real(dp), intent(in) :: c
     real(dp), intent(in) :: phi
     intgrandFast=exp(-c/cos(phi))
     return
  end function

END MODULE mo_neutrons
