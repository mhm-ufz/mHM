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
  !>        \param[in] "real(dp), dimension(:,:) :: sm" Soil Moisture
  !>        \param[in] "real(dp), dimension(:)   :: Horizons" Horizon depths
  !>        \param[in] "real(dp)                 :: N0" dry neutron counts
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(size(sm,1)) :: neutrons" Neutron counts
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
  !         N0=1500cph, sm(1,1)=700mm, Horizons(1)=200mm
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

  subroutine DesiletsN0(sm, Horizons, N0, neutrons)

    use mo_mhm_constants, only: Desilets_a0, Desilets_a1, Desilets_a2
    implicit none
    
    real(dp), dimension(:,:),        intent(in)  :: sm
    real(dp), dimension(:),          intent(in)  :: Horizons
    real(dp),                        intent(in)  :: N0          ! from global parameters
    real(dp), dimension(size(sm,1)), intent(out) :: neutrons
    
    ! only use first soil layer
    neutrons(:) = N0 * ( Desilets_a1 + Desilets_a0 / (sm(:,1)/Horizons(1) + Desilets_a2))
  
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
  !                          COSMIC-parameterset, output(cells) )
  !
  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:,:) :: sm" Soil Moisture
  !>        \param[in] "real(dp), dimension(:)   :: Horizons" Horizon depths
  !>        \param[in] "real(dp), dimension(:)   :: params" ! N0, N1, N2, alpha0, alpha1, L30, L31
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(size(sm,1)) :: neutrons" Neutron counts
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
  
  subroutine COSMIC(sm, Horizons, params, neutrons)
    
    use mo_mhm_constants, only: H2Odens, &
        COSMIC_bd, COSMIC_vwclat, COSMIC_N, COSMIC_alpha, &
        COSMIC_L1, COSMIC_L2, COSMIC_L3, COSMIC_L4
    use mo_constants, only: PI_dp
    implicit none
    
    real(dp), dimension(:,:),        intent(in)  :: sm
    real(dp), dimension(:),          intent(in)  :: Horizons
    real(dp), dimension(:),          intent(in)  :: params ! 1: N0, 2: N1, 3: N2, 4: alpha0, 5: alpha1, 6: L30, 7. L31
    real(dp), dimension(size(sm,1)), intent(out) :: neutrons

    ! local variables
    real(dp) :: zdeg         
    real(dp) :: zrad         
    real(dp) :: ideg         
    real(dp) :: costheta     
    real(dp) :: dtheta       
   
    real(dp), dimension(size(Horizons))   :: dz          ! Soil layers (cm)
    real(dp), dimension(size(Horizons))   :: zthick      ! Soil layer thickness (cm)
    real(dp), dimension(:,:), allocatable :: wetsoidens  ! Density of wet soil layer (g/cm3)
    real(dp), dimension(:,:), allocatable :: wetsoimass  ! Mass of wet soil layer (g)
    real(dp), dimension(:,:), allocatable :: isoimass    ! Integrated dry soil mass above layer (g)
    real(dp), dimension(:,:), allocatable :: iwatmass    ! Integrated water mass above layer (g)
    real(dp), dimension(:,:), allocatable :: iwetsoimass ! Integrated wet soil mass above layer (g)
    real(dp), dimension(:,:), allocatable :: hiflux      ! High energy neutron flux
    real(dp), dimension(:,:), allocatable :: fastpot     ! Fast neutron source strength of layer
    real(dp), dimension(:,:), allocatable :: h2oeffdens  ! "Effective" density of water in layer (g/cm3)
    real(dp), dimension(:,:), allocatable :: h2oeffmass  ! "Effective" mass of water in layer (g)
    real(dp), dimension(:,:), allocatable :: ih2oeffmass ! Integrated water mass above layer (g)
    real(dp), dimension(:,:), allocatable :: idegrad     ! Integrated neutron degradation factor (-)
    real(dp), dimension(:,:), allocatable :: fastflux    ! Contribution to above-ground neutron flux
    real(dp), dimension(:)  , allocatable :: totflux     ! Total flux of above-ground fast neutrons
    real(dp), dimension(:,:), allocatable :: normfast    ! Normalized contribution to neutron flux (-) [weighting factors]
    real(dp), dimension(:,:), allocatable :: inormfast   ! Cumulative fraction of neutrons (-)
    !
    integer(i4) :: layers=1                 ! Total number of soil layers
    integer(i4) :: profiles=1               ! Total number of soil moisture profiles
    integer(i4) :: ll=1,pp=0
    integer(i4) :: angle, angledz, maxangle ! loop indices for an integration interval
    !
    layers   = size(sm,2) ! 2
    profiles = size(sm,1) ! 34
    
    allocate(totflux(profiles))
    allocate(wetsoidens(layers,profiles),wetsoimass(layers,profiles),&
             iwetsoimass(layers,profiles),hiflux(layers,profiles),fastpot(layers,profiles),&
             h2oeffdens(layers,profiles),h2oeffmass(layers,profiles),ih2oeffmass(layers,profiles),&
             idegrad(layers,profiles),fastflux(layers,profiles),normfast(layers,profiles),&
             inormfast(layers,profiles),isoimass(layers,profiles),iwatmass(layers,profiles))

    dz(:)            = 0.0_dp * params(1) ! <-- this multiplication with params(1) is not needed, only to make params USED
    !                                     !     PLEASE remove when possible 
    zthick(:)        = 0.0_dp
    wetsoidens(:,:)  = 0.0_dp
    wetsoimass(:,:)  = 0.0_dp
    iwetsoimass(:,:) = 0.0_dp
    isoimass(:,:)    = 0.0_dp
    iwatmass(:,:)    = 0.0_dp
    hiflux(:,:)      = 0.0_dp
    fastpot(:,:)     = 0.0_dp
    h2oeffdens(:,:)  = 0.0_dp
    h2oeffmass(:,:)  = 0.0_dp
    ih2oeffmass(:,:) = 0.0_dp
    idegrad(:,:)     = 0.0_dp
    fastflux(:,:)    = 0.0_dp
    normfast(:,:)    = 0.0_dp
    inormfast(:,:)   = 0.0_dp 
    totflux(:)       = 0.0_dp
    
    ! Soil Layers and Thicknesses are constant in mHM, they could be defined outside of this function
    dz(:) = Horizons(:)/10.0_dp ! from mm to cm
    zthick(1) = dz(1) - 0.0_dp
    do ll = 2,layers
       zthick(ll) = dz(ll) - dz(ll-1)
    enddo
    
    !
    ! Angle distribution parameters (HARDWIRED)
    ! rr: Using 0.5 deg angle intervals appears to be sufficient
    ! rr: (smaller angles increase the computing time for COSMIC)
    ideg     = 0.5_dp                   ! ideg ultimately controls the number of trips through
    angledz  = nint(ideg*10.0_dp,i4)    ! the ANGLE loop. Make sure the 10.0 is enough
    maxangle = 900_i4 - angledz         ! to create integers with no remainder
    dtheta   = ideg*(PI_dp/180.0_dp)

    !
    do pp = 1,profiles
       do ll = 1,layers
          
          ! High energy neutron downward flux
          ! The integration is now performed at the node of each layer (i.e., center of the layer)
          h2oeffdens(ll,pp) = ((sm(pp,ll) / zthick(ll) / 10.0_dp +COSMIC_vwclat)*H2Odens)/1000.0_dp  

          ! Assuming an area of 1 cm2
          if(ll > 1) then
             isoimass(ll,pp) = isoimass(ll-1,pp) + COSMIC_bd*(0.5_dp*zthick(ll-1))*1.0_dp &
                                             + COSMIC_bd*(0.5_dp*zthick(ll))*1.0_dp
             iwatmass(ll,pp) = iwatmass(ll-1,pp) + h2oeffdens(ll-1,pp)*(0.5_dp*zthick(ll-1))*1.0_dp &
                                             + h2oeffdens(ll,pp)*(0.5_dp*zthick(ll))*1.0_dp
          else
             isoimass(ll,pp) = COSMIC_bd*(0.5_dp*zthick(ll))*1.0_dp 
             iwatmass(ll,pp) = h2oeffdens(ll,pp)*(0.5_dp*zthick(ll))*1.0_dp
          endif

          hiflux(ll,pp)  = COSMIC_N*exp(-(isoimass(ll,pp)/COSMIC_L1 + iwatmass(ll,pp)/COSMIC_L2) )
          fastpot(ll,pp) = zthick(ll)*hiflux(ll,pp)*(COSMIC_alpha*COSMIC_bd + h2oeffdens(ll,pp))

          ! This second loop needs to be done for the distribution of angles for fast neutron release
          ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
          ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.  

          do angle=0,maxangle,angledz
             zdeg     = real(angle,dp)/10.0_dp   ! 0.0  0.5  1.0  1.5 ...
             zrad     = (zdeg*PI_dp)/180.0_dp
             costheta = cos(zrad)

             ! Angle-dependent low energy (fast) neutron upward flux
             fastflux(ll,pp) = fastflux(ll,pp) + fastpot(ll,pp) * &
                               exp(-(isoimass(ll,pp)/COSMIC_L3 + iwatmass(ll,pp)/COSMIC_L4)/costheta)*dtheta
          enddo

          ! After contribution from all directions are taken into account,
          ! need to multiply fastflux by 2/pi
          fastflux(ll,pp) = (2.0_dp/PI_dp)*fastflux(ll,pp)

          ! Low energy (fast) neutron upward flux
          totflux(pp) = totflux(pp) + fastflux(ll,pp)

       enddo
    enddo
    
    neutrons = totflux(:)

    deallocate(totflux, wetsoidens, wetsoimass, iwetsoimass, hiflux,&
           fastpot, h2oeffdens, h2oeffmass, ih2oeffmass, idegrad, fastflux,&
           normfast, inormfast, isoimass, iwatmass)
           
  end subroutine COSMIC
  
END MODULE mo_neutrons
