!> \file mo_neutrons.f90
!> \brief Models to predict neutron intensities above soils
!> \authors Martin Schroen
!> \date Mar 2015

! THIS MODULE IS WORK IN PROGRESS, DO NOT USE FOR RESEARCH.

! TODO make it faster with pre-calculated horizons and variables
! TODO outsource constants to mo_constants
! TODO use global parameters as linear model 

MODULE mo_neutrons

  USE mo_kind, ONLY: i4, sp, dp
  IMPLICIT NONE

  PUBLIC :: COSMIC     ! Neutron forward model using particle transport physics, see Shuttleworth et al. 2013
  PUBLIC :: DesiletsN0 ! inverse \theta(N) relation based on Desilets et al. 2010

  PRIVATE
  !

CONTAINS
  
! ------------------------------------------------------------------
! N0 formula based on Desilets et al. 2010                          
! ------------------------------------------------------------------
  
  subroutine DesiletsN0(sm, Horizons, N0, neutrons)
  
    use mo_mhm_constants, only: a0, a1, a2
	implicit none
	real(dp), dimension(:,:), intent(in)       :: sm
	real(dp), dimension(:), intent(in)         :: Horizons
	real(dp), intent(in)                       :: N0 ! from global parameters
	real(dp), dimension(size(sm,1)), intent(out) :: neutrons
	
	! only use first soil layer
	neutrons(:) = N0 * ( a1 + a0 / (sm(:,1)/Horizons(1) + a2))
  
  end subroutine DesiletsN0

! ------------------------------------------------------------------
! COSMIC model based on Shuttleworth et al. 2013                    
! ------------------------------------------------------------------

  subroutine COSMIC(sm, Horizons, params, neutrons)
  
    use mo_mhm_constants, only: H2Odens
    implicit none
	real(dp), dimension(:,:), intent(in)         :: sm
	real(dp), dimension(:), intent(in)           :: Horizons
	real(dp), dimension(:), intent(in)           :: params ! 1: N0, 2: N1, 3: N2, 4: alpha0, 5: alpha1, 6: L30, 7. L31
	real(dp), dimension(size(sm,1)), intent(out) :: neutrons
	!
    real(dp)            :: bd=1.4020_dp          ! Dry soil bulk density (g/m3)
	real(dp)            :: vwclat=0.0753_dp       ! Volumetric "lattice" water content (m3/m3)
	real(dp)            :: N=510.5173790200_dp            ! High energy neutron flux (-)
	real(dp)            :: alpha=0.2392421548_dp        ! Ratio of Fast Neutron Creation Factor (Soil to Water), alpha (-)
	real(dp)            :: L1=161.98621864_dp           ! High Energy Soil Attenuation Length (g/cm2)
	real(dp)            :: L2=129.14558985_dp           ! High Energy Water Attenuation Length (g/cm2)
	real(dp)            :: L3=107.8220456200_dp           ! Fast Neutron Soil Attenuation Length (g/cm2)
	real(dp)            :: L4=3.1627190566_dp           ! Fast Neutron Water Attenuation Length (g/cm2)
	real(dp)            :: zdeg         
	real(dp)            :: zrad         
	real(dp)            :: ideg         
	real(dp)            :: costheta     
	real(dp)            :: dtheta       
	real(dp), parameter :: pi=3.14159265359_dp
  
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
	integer(i4)            :: nlayers=1                 ! Total number of soil layers
	integer(i4)            :: nprof=1                   ! Total number of soil moisture profiles
	integer(i4)            :: i=1,j=0
	integer(i4)            :: angle, angledz, maxangle  ! loop indices for an integration interval
	!
	nlayers = size(sm,2) ! 2
	nprof   = size(sm,1) ! 34
	
	allocate(totflux(nprof))
	allocate(wetsoidens(nlayers,nprof),wetsoimass(nlayers,nprof),&
			 iwetsoimass(nlayers,nprof),hiflux(nlayers,nprof),fastpot(nlayers,nprof),&
			 h2oeffdens(nlayers,nprof),h2oeffmass(nlayers,nprof),ih2oeffmass(nlayers,nprof),&
			 idegrad(nlayers,nprof),fastflux(nlayers,nprof),normfast(nlayers,nprof),&
			 inormfast(nlayers,nprof),isoimass(nlayers,nprof),iwatmass(nlayers,nprof))

	do i = 1,nlayers
	   dz(i)     = 0.0_dp
	   zthick(i) = 0.0_dp
	   do j = 1,nprof
		  wetsoidens(i,j)  = 0.0_dp
		  wetsoimass(i,j)  = 0.0_dp
		  iwetsoimass(i,j) = 0.0_dp
		  isoimass(i,j)    = 0.0_dp
		  iwatmass(i,j)    = 0.0_dp
		  hiflux(i,j)      = 0.0_dp
		  fastpot(i,j)     = 0.0_dp
		  h2oeffdens(i,j)  = 0.0_dp
		  h2oeffmass(i,j)  = 0.0_dp
		  ih2oeffmass(i,j) = 0.0_dp
		  idegrad(i,j)     = 0.0_dp
		  fastflux(i,j)    = 0.0_dp
		  normfast(i,j)    = 0.0_dp
		  inormfast(i,j)   = 0.0_dp 
		  !if(nlayers == 1)
		  totflux(j) = 0.0_dp
	   enddo
	enddo
	
	! Soil Layers and Thicknesses are constant in mHM, they could be defined outside of this function
	dz(:) = Horizons(:)/10.0_dp ! from mm to cm
	zthick(1) = dz(1) - 0.0_dp
	do i = 2,nlayers
	   zthick(i) = dz(i) - dz(i-1)
	enddo
	
	!
	! Angle distribution parameters (HARDWIRED)
	!rr: Using 0.5 deg angle intervals appears to be sufficient
	!rr: (smaller angles increase the computing time for COSMIC)
	ideg     = 0.5_dp                   ! ideg ultimately controls the number of trips through
	angledz  = nint(ideg*10.0_dp)       ! the ANGLE loop. Make sure the 10.0 is enough
	maxangle = 900_dp - angledz         ! to create integers with no remainder
	dtheta   = ideg*(PI/180.0_dp)

	if ( real(angledz) /= ideg*10.0_dp ) &
		write(*,*) 'ideg*10.0 must result in an integer - it results in ',ideg*10.0_dp
	!
	do j = 1,nprof
	   do i = 1,nlayers
		  
		  ! High energy neutron downward flux
		  ! The integration is now performed at the node of each layer (i.e., center of the layer)
		  h2oeffdens(i,j) = ((sm(j,i) / zthick(i) / 10.0_dp +vwclat)*H2Odens)/1000.0_dp  

		  ! Assuming an area of 1 cm2
		  if(i > 1) then
			 isoimass(i,j) = isoimass(i-1,j) + bd*(0.5_dp*zthick(i-1))*1.0_dp &
											 + bd*(0.5_dp*zthick(i  ))*1.0_dp
			 iwatmass(i,j) = iwatmass(i-1,j) + h2oeffdens(i-1,j)*(0.5_dp*zthick(i-1))*1.0_dp &
											 + h2oeffdens(i  ,j)*(0.5_dp*zthick(i  ))*1.0_dp
		  else
			 isoimass(i,j) = bd*(0.5_dp*zthick(i))*1.0_dp 
			 iwatmass(i,j) = h2oeffdens(i,j)*(0.5_dp*zthick(i))*1.0_dp
		  endif

		  hiflux(i,j)  = N*exp(-(isoimass(i,j)/L1 + iwatmass(i,j)/L2) )
		  fastpot(i,j) = zthick(i)*hiflux(i,j)*(alpha*bd + h2oeffdens(i,j))

		  ! This second loop needs to be done for the distribution of angles for fast neutron release
		  ! the intent is to loop from 0 to 89.5 by 0.5 degrees - or similar.
		  ! Because Fortran loop indices are integers, we have to divide the indices by 10 - you get the idea.  

		  do angle=0,maxangle,angledz
			 zdeg     = real(angle)/10.0_dp   ! 0.0  0.5  1.0  1.5 ...
			 zrad     = (zdeg*pi)/180.0_dp
			 costheta = cos(zrad)

			 ! Angle-dependent low energy (fast) neutron upward flux
			 fastflux(i,j) = fastflux(i,j) + fastpot(i,j)*exp(-(isoimass(i,j)/L3 + iwatmass(i,j)/L4)/costheta)*dtheta
		  enddo

		  ! After contribution from all directions are taken into account,
		  ! need to multiply fastflux by 2/PI
		  fastflux(i,j) = (2.0_dp/pi)*fastflux(i,j)

		  ! Low energy (fast) neutron upward flux
		  totflux(j) = totflux(j) + fastflux(i,j)

	   enddo
!	   if (.not. totflux(j) .gt. 0.01_dp ) write(*,*) totflux(j)
	enddo
	
	neutrons = totflux(:)

	deallocate(totflux, wetsoidens, wetsoimass, iwetsoimass, hiflux,&
           fastpot, h2oeffdens, h2oeffmass, ih2oeffmass, idegrad, fastflux,&
           normfast, inormfast, isoimass, iwatmass)
		   
  end subroutine COSMIC
  
END MODULE mo_neutrons
