!>       \file mo_mrm_mpr.f90

!>       \brief Perform Multiscale Parameter Regionalization on Routing Parameters

!>       \details This module contains the subroutine for calculating the regionalized
!>       routing parameters (beta-parameters) given the five global routing parameters
!>       (gamma) at the level 0 scale.

!>       \authors Luis Samaniego, Stephan Thober

!>       \date Aug 2015

! Modifications:

module mo_mrm_mpr
  use mo_kind, only : dp
  implicit none
  public :: reg_rout
  private
contains

  ! ----------------------------------------------------------------------------

  !    NAME
  !        reg_rout

  !    PURPOSE
  !>       \brief Regionalized routing

  !>       \details sets up the Regionalized Routing parameters
  !>       Global parameters needed (see mhm_parameter.nml):
  !>       - param(1) = muskingumTravelTime_constant    
  !>       - param(2) = muskingumTravelTime_riverLength 
  !>       - param(3) = muskingumTravelTime_riverSlope  
  !>       - param(4) = muskingumTravelTime_impervious  
  !>       - param(5) = muskingumAttenuation_riverSlope 

  !    INTENT(IN)
  !>       \param[in] "real(dp), dimension(5) :: param"  input parameter
  !>       \param[in] "real(dp), dimension(:) :: length" [m] total length
  !>       \param[in] "real(dp), dimension(:) :: slope"  average slope
  !>       \param[in] "real(dp), dimension(:) :: fFPimp" fraction of the flood plain withimpervious layer
  !>       \param[in] "real(dp) :: TS"                   - [h] time step in

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: C1" routing parameter C1 (Chow, 25-41)
  !>       \param[out] "real(dp), dimension(:) :: C2" routing parameter C2 (")

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:

  subroutine reg_rout(param, length, slope, fFPimp, TS, C1, C2)
    implicit none

    ! input parameter
    real(dp), dimension(5), intent(in) :: param

    ! [m] total length
    real(dp), dimension(:), intent(in) :: length

    ! average slope
    real(dp), dimension(:), intent(in) :: slope

    ! fraction of the flood plain withimpervious layer
    real(dp), dimension(:), intent(in) :: fFPimp

    ! - [h] time step in
    real(dp), intent(in) :: TS

    ! routing parameter C1 (Chow, 25-41)
    real(dp), dimension(:), intent(out) :: C1

    ! routing parameter C2 (")
    real(dp), dimension(:), intent(out) :: C2

    ! stream slope max
    real(dp) :: ssMax

    ! [d] Muskingum travel time parameter
    real(dp), dimension(size(fFPimp, 1)) :: K

    ! [1] Muskingum diffusion parameter (attenuation)
    real(dp), dimension(size(fFPimp, 1)) :: xi


    ! normalize stream bed slope
    ssMax = maxval(slope(:))

    ! New regional relationship; K = f(length, slope, & fFPimp)
    K = param(1) + param(2) * (length * 0.001_dp) &
            + param(3) * slope &
            + param(4) * fFPimp

    ! Xi = f(slope)
    xi = param(5) * (1.0_dp + slope / ssMax)

    ! constraints on Xi
    xi = merge(0.5_dp, xi, xi > 0.5_dp)
    xi = merge(0.005_dp, xi, xi < 0.005_dp)

    ! constrains on Ki
    K = merge(0.5_dp * TS / xi, K, K > 0.5_dp * TS / xi)
    K = merge(0.5_dp * TS / (1.0_dp - xi), K, K < 0.5_dp * TS / (1.0_dp - xi))

    ! Muskingum parameters
    C1 = TS / (K * (1.0_dp - xi) + 0.5_dp * TS)
    C2 = 1.0_dp - C1 * K / TS

  end subroutine reg_rout
end module mo_mrm_mpr
