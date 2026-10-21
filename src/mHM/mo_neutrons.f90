!> \file mo_neutrons.f90
!> \copydoc mo_neutrons

!> \brief Models to predict neutron intensities above soils
!> \details The number of neutrons above the ground is directly related to
!! the number soil water content in the ground, air, vegetation and/or snow.
!! This module forward-models neutron abundance as a state variable for each cell.
!> \authors Martin Schroen
!> \date Mar 2015
!> \warning THIS MODULE IS WORK IN PROGRESS, DO NOT USE FOR RESEARCH.
!> \version 0.1
!> \changelog
!! - 2026: component ownership, explicit interfaces, and retained-kernel documentation cleanup.
!! - 2026: read-only geometry and scalar accumulators replace per-thread neutron work arrays.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_neutrons

  ! One immutable geometry cache is shared by the cell calculations of each container.

  USE mo_kind, ONLY: i4, dp
  IMPLICIT NONE

  ! Neutron forward model using particle transport physics, see Shuttleworth et al. 2013
  PUBLIC :: COSMIC

  ! inverse \theta(N) relation based on Desilets et al. 2010
  PUBLIC :: DesiletsN0

  ! integration tabular for approximating the neutron flux integral
  PUBLIC :: TabularIntegralAFast

  PRIVATE

  !> \class neutron_geometry_t
  !> \brief Initialization-time geometry shared read-only by neutron calculations.
  type, public :: neutron_geometry_t
    private
    real(dp), allocatable :: thickness(:) !< Simulated soil-layer thicknesses [mm].
    integer(i4), allocatable :: first_sample(:) !< First 10 mm sample of each layer.
    integer(i4), allocatable :: last_sample(:) !< Last 10 mm sample of each layer.
  contains
    procedure :: prepare => neutron_geometry_prepare
    procedure :: destroy => neutron_geometry_destroy
  end type neutron_geometry_t

CONTAINS

  !> \brief Cache invariant layer thicknesses and sampling ranges during initialization.
  subroutine neutron_geometry_prepare(self, horizons)
    class(neutron_geometry_t), intent(inout), target :: self
    real(dp), intent(in) :: horizons(:) !< Simulated horizon bottoms [mm], excluding the last soil layer.
    integer(i4) :: layer, intervals
    real(dp) :: lower
    call self%destroy()
    allocate(self%thickness(size(horizons)), self%first_sample(size(horizons)), self%last_sample(size(horizons)))
    intervals = nint(maxval(horizons)/10.0_dp)
    lower = 0.0_dp
    do layer = 1, size(horizons)
      self%thickness(layer) = horizons(layer) - lower
      call neutron_sample_range(lower, horizons(layer), intervals, self%first_sample(layer), self%last_sample(layer))
      lower = horizons(layer)
    end do
  end subroutine neutron_geometry_prepare

  !> \brief Release the owning container's cached geometry.
  subroutine neutron_geometry_destroy(self)
    class(neutron_geometry_t), intent(inout), target :: self
    if (allocated(self%thickness)) deallocate(self%thickness)
    if (allocated(self%first_sample)) deallocate(self%first_sample)
    if (allocated(self%last_sample)) deallocate(self%last_sample)
  end subroutine neutron_geometry_destroy

  !> \brief Locate samples using the original half-open layer intervals and 10 mm spacing.
  pure subroutine neutron_sample_range(lower, upper, intervals, first, last)
    real(dp), intent(in) :: lower !< Layer top [mm].
    real(dp), intent(in) :: upper !< Layer bottom [mm].
    integer(i4), intent(in) :: intervals !< Number of samples over the simulated profile.
    integer(i4), intent(out) :: first !< First sample index; one for an empty range.
    integer(i4), intent(out) :: last !< Last sample index; zero for an empty range.
    integer(i4) :: sample
    real(dp) :: depth
    first = 1_i4
    last = 0_i4
    depth = 0.0_dp
    do sample = 1, intervals
      if (depth >= lower .and. depth < upper) then
        if (last == 0_i4) first = sample
        last = sample
      end if
      depth = depth + 10.0_dp
    end do
  end subroutine neutron_sample_range

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

  recursive subroutine DesiletsN0(SoilMoisture, Horizon_depth, Bd, latWater, N0, neutrons, geometry)
    use mo_mhm_constants, only: Desilets_a0, Desilets_a1, Desilets_a2
    implicit none
    real(dp), intent(in) :: SoilMoisture(:) !< Soil water storage in the simulated horizons [mm].
    real(dp), intent(in) :: Horizon_depth(:) !< Bottom depths of the simulated horizons [mm].
    real(dp), intent(in) :: Bd(:) !< Bulk density [g cm-3].
    real(dp), intent(in) :: latWater(:) !< Lattice water storage [mm].
    real(dp), intent(in) :: N0 !< Calibrated dry neutron count.
    real(dp), intent(inout) :: neutrons !< Simulated neutron count.
    type(neutron_geometry_t), optional, intent(in) :: geometry !< Shared initialization-time geometry.

    integer(i4) :: nLayers, LL, nn, first, last
    real(dp) :: average_swc, average_bd, D_86_in_cm, D_86_in_mm
    real(dp) :: lower, thickness, depth, weight_10mm_spacing, grav_swc, layer_weight, total_weight

    nLayers = size(SoilMoisture)
    average_swc = 0.0_dp
    lower = 0.0_dp
    do LL = 1, nLayers
      thickness = Horizon_depth(LL) - lower
      if (present(geometry)) thickness = geometry%thickness(LL)
      average_swc = average_swc + SoilMoisture(LL) / thickness
      lower = Horizon_depth(LL)
    end do
    average_swc = average_swc / real(nLayers, dp)
    average_bd = sum(Bd) / real(size(Bd), dp)

    ! Moisture-dependent penetration depth [cm and mm], recomputed for each cell and interval.
    D_86_in_cm = (1.0_dp/average_bd) * &
      (8.321_dp + 0.14249_dp * (0.96655_dp + exp(-0.01_dp)) * (20.0_dp + average_swc)/(0.0429_dp + average_swc))
    D_86_in_mm = D_86_in_cm * 10.0_dp

    ! Accumulate each layer's samples in their original order, then combine layers in order.
    average_swc = 0.0_dp
    total_weight = 0.0_dp
    lower = 0.0_dp
    do LL = 1, nLayers
      thickness = Horizon_depth(LL) - lower
      if (present(geometry)) then
        thickness = geometry%thickness(LL)
        first = geometry%first_sample(LL)
        last = geometry%last_sample(LL)
      else
        call neutron_sample_range(lower, Horizon_depth(LL), nint(maxval(Horizon_depth)/10.0_dp), first, last)
      end if
      layer_weight = 0.0_dp
      depth = real(first - 1_i4, dp) * 10.0_dp
      do nn = first, last
        weight_10mm_spacing = exp(-2.0_dp * depth / D_86_in_mm)
        layer_weight = layer_weight + weight_10mm_spacing
        depth = depth + 10.0_dp
      end do
      grav_swc = (SoilMoisture(LL) + latWater(LL))/thickness/Bd(LL)
      average_swc = average_swc + grav_swc * layer_weight
      total_weight = total_weight + layer_weight
      lower = Horizon_depth(LL)
    end do
    average_swc = average_swc / total_weight
    neutrons = N0 * (Desilets_a1 + Desilets_a0 / (average_swc + Desilets_a2))
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
  recursive subroutine COSMIC(SoilMoisture, Horizons, neutron_integral_AFast, &
      interc, snowpack, L1_N0, L1_bulkDens, L1_latticeWater, L1_COSMICL3, neutrons, geometry)
    use mo_mhm_constants, only: H2Odens, COSMIC_alpha, COSMIC_L1, COSMIC_L2, COSMIC_L4
    use mo_constants, only: PI_dp
    implicit none
    real(dp), intent(in) :: SoilMoisture(:) !< Soil water storage in the simulated horizons [mm].
    real(dp), intent(in) :: Horizons(:) !< Bottom depths of the simulated horizons [mm].
    real(dp), intent(in) :: neutron_integral_AFast(:) !< Cached angular neutron-flux integral table.
    real(dp), intent(in) :: interc !< Canopy interception storage [mm].
    real(dp), intent(in) :: snowpack !< Snow water storage [mm].
    real(dp), intent(in) :: L1_N0 !< Calibrated dry neutron count.
    real(dp), intent(in) :: L1_bulkDens(:) !< Horizon bulk density [g cm-3].
    real(dp), intent(in) :: L1_latticeWater(:) !< Horizon lattice water storage [mm].
    real(dp), intent(in) :: L1_COSMICL3(:) !< Horizon COSMIC attenuation length [g cm-2].
    real(dp), intent(inout) :: neutrons !< Simulated neutron count.
    type(neutron_geometry_t), optional, intent(in) :: geometry !< Shared initialization-time geometry.

    integer(i4) :: layer
    real(dp) :: lower, thickness, previous_thickness, previous_soil_mass, previous_water_mass, previous_water_density
    real(dp) :: soil_mass, water_mass, water_height, water_density, bd, lambdaHigh, lambdaFast, hiflux, xeff, fastflux, totflux

    ! Preserve the existing disabled surface-layer contribution, including its thickness in the next soil mass.
    previous_thickness = (snowpack + interc) / 10.0_dp
    previous_soil_mass = 0.0_dp
    previous_water_mass = 0.0_dp
    previous_water_density = 0.0_dp
    lower = 0.0_dp
    totflux = 0.0_dp
    do layer = 1, size(SoilMoisture)
      thickness = Horizons(layer) - lower
      if (present(geometry)) thickness = geometry%thickness(layer)
      thickness = thickness / 10.0_dp
      soil_mass = 0.0_dp
      water_mass = 0.0_dp
      water_density = 0.0_dp
      if (thickness > 0.0_dp) then
        bd = L1_bulkDens(layer)
        water_height = SoilMoisture(layer) / 10.0_dp
        water_density = (water_height + L1_latticeWater(layer)/10.0_dp)/thickness*H2Odens/1000.0_dp

        ! Midpoint column masses depend only on the preceding layer; retain the original operation order.
        soil_mass = bd*(0.5_dp*thickness)*1.0_dp
        water_mass = water_density*(0.5_dp*thickness)*1.0_dp
        soil_mass = soil_mass + previous_soil_mass + bd*(0.5_dp*previous_thickness)*1.0_dp
        water_mass = water_mass + previous_water_mass + previous_water_density*(0.5_dp*previous_thickness)*1.0_dp
        lambdaHigh = soil_mass/COSMIC_L1 + water_mass/COSMIC_L2
        lambdaFast = soil_mass/L1_COSMICL3(layer) + water_mass/COSMIC_L4
        hiflux = exp(-lambdaHigh)
        xeff = thickness*(COSMIC_alpha*bd + water_density)
        call lookUpIntegral(fastflux, neutron_integral_AFast, lambdaFast)
        fastflux = (2.0_dp/PI_dp) * fastflux
        totflux = totflux + hiflux * xeff * fastflux
      end if
      previous_thickness = thickness
      previous_soil_mass = soil_mass
      previous_water_mass = water_mass
      previous_water_density = water_density
      lower = Horizons(layer)
    end do
    neutrons = L1_N0*totflux
  end subroutine COSMIC

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
     real(dp)                                     :: res !< Computed integral or accumulated quadrature result.
     real(dp), external                           :: f !< Integrand evaluated as f(c, x).
     real(dp),                        intent(in)  :: c !< Dimensionless neutron attenuation parameter.
     real(dp),                        intent(in)  :: xmax !< Upper integration limit.
     real(dp),                        intent(in)  :: xmin !< Lower integration limit.
     real(dp),                        optional    :: eps !< Requested quadrature error tolerance.
     integer(i4),                     optional    :: steps !< Maximum number of quadrature subdivisions.
     real(dp),                        optional    :: fxmin !< Integrand value at the lower limit; computed when absent.
     real(dp),                        optional    :: fxmax !< Integrand value at the upper limit; computed when absent.

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
     real(dp)                                     :: res !< Computed integral or accumulated quadrature result.
     real(dp), external                           :: f !< Integrand evaluated as f(c, x).
     real(dp),                        intent(in)  :: c !< Dimensionless neutron attenuation parameter.
     real(dp),                        intent(in)  :: xmax !< Upper integration limit.
     real(dp),                        intent(in)  :: xmin !< Lower integration limit.
     real(dp),                        intent(in)  :: eps !< Requested quadrature error tolerance.
     integer(i4),                     intent(in)  :: steps !< Maximum number of quadrature subdivisions.
     real(dp),                        intent(in)  :: fxmin !< Integrand value at the lower integration limit.
     real(dp),                        intent(in)  :: fxmax !< Integrand value at the upper integration limit.

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
     real(dp)                                     :: res !< Computed integral or accumulated quadrature result.
     real(dp), external                           :: f !< Integrand evaluated as f(c, x).
     real(dp),                        intent(in)  :: c !< Dimensionless neutron attenuation parameter.
     real(dp),                        intent(in)  :: xmax !< Upper integration limit.
     real(dp),                        intent(in)  :: xmin !< Lower integration limit.
     real(dp),                        intent(in)  :: eps !< Requested quadrature error tolerance.
     real(dp),                        intent(in)  :: fxmin !< Integrand value at the lower integration limit.
     real(dp),                        intent(in)  :: fxmax !< Integrand value at the upper integration limit.

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
     real(dp), dimension(:)              :: integral !< Tabulated angular integral and its lookup metadata.
     real(dp), intent(in)                :: maxC !< Maximum tabulated attenuation parameter.

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
     real(dp), dimension(:)              :: integral !< Tabulated angular integrals and scaled derivatives.
     integer(i4), intent(in)             :: intsize !< Number of tabulated interpolation intervals.
     real(dp), intent(in)                :: maxC !< Maximum tabulated attenuation parameter.

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
  recursive subroutine lookUpIntegral(res,integral,c)
     use mo_constants, only: PI_dp
     implicit none
     real(dp),              intent(out):: res !< Computed integral or accumulated quadrature result.
     real(dp), dimension(:),intent(in) :: integral !< Tabulated angular integral and its lookup metadata.
     real(dp),              intent(in) :: c !< Dimensionless neutron attenuation parameter.

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
     real(dp)                         :: res !< Computed integral or accumulated quadrature result.
     real(dp), dimension(:),intent(in):: integral !< Tabulated angular integral and its lookup metadata.
     integer(i4), intent(in)          :: intsize !< Number of tabulated interpolation intervals.
     real(dp), intent(in)             :: c !< Dimensionless neutron attenuation parameter.
     real(dp), intent(in)             :: maxC !< Maximum tabulated attenuation parameter.

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
     real(dp)                                     :: res !< Computed integral or accumulated quadrature result.
     real(dp),                        intent(in)  :: x !< Dimensionless evaluation point.
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
     real(dp)                                     :: res !< Computed integral or accumulated quadrature result.
     real(dp),                        intent(in)  :: c !< Dimensionless neutron attenuation parameter.

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
     real(dp), intent(in) :: c !< Dimensionless neutron attenuation parameter.
     real(dp), intent(in) :: phi !< Zenith angle [radians].
     intgrandFast=exp(-c/cos(phi))
     return
  end function

  function intgrandDerivFast(c,phi)
     implicit none
     real(dp) :: intgrandDerivFast
     real(dp), intent(in) :: c !< Dimensionless neutron attenuation parameter.
     real(dp), intent(in) :: phi !< Zenith angle [radians].
     intgrandDerivFast=(-1.0_dp/cos(phi))*exp(-c/cos(phi))
     return
  end function

  function expPolynomDeg3(x,a,b,c,d)
     implicit none
     real(dp) :: expPolynomDeg3
     real(dp), intent(in) :: x !< Dimensionless evaluation point.
     real(dp), intent(in) :: a !< Cubic coefficient in the exponential polynomial.
     real(dp), intent(in) :: b !< Quadratic coefficient in the exponential polynomial.
     real(dp), intent(in) :: c !< Linear coefficient in the exponential polynomial.
     real(dp), intent(in) :: d !< Constant coefficient in the exponential polynomial.

     expPolynomDeg3=exp(a*x**3+b*x**2+c*x+d)
     return
  end function

  ! hermite polynoms
  function h00(t)
     implicit none
     real(dp)             :: h00
     real(dp), intent(in) :: t !< Normalized Hermite interpolation coordinate.
     h00=2.0_dp*t**3-3.0_dp*t**2+1.0_dp
     return
  end function

  function h01(t)
     implicit none
     real(dp)             :: h01
     real(dp), intent(in) :: t !< Normalized Hermite interpolation coordinate.
     h01=-2.0_dp*t**3+3.0_dp*t**2
     return
  end function

  function h10(t)
     implicit none
     real(dp)             :: h10
     real(dp), intent(in) :: t !< Normalized Hermite interpolation coordinate.
     h10=t**3-2.0_dp*t**2+t
     return
  end function

  function h11(t)
     implicit none
     real(dp)             :: h11
     real(dp), intent(in) :: t !< Normalized Hermite interpolation coordinate.
     h11=t**3-t**2
     return
  end function

END MODULE mo_neutrons
