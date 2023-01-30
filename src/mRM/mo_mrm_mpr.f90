!> \file mo_mrm_mpr.f90
!> \brief \copybrief mo_mrm_mpr
!> \details \copydetails mo_mrm_mpr

!> \brief Perform Multiscale Parameter Regionalization on Routing Parameters
!> \details This module contains the subroutine for calculating the regionalized
!! routing parameters (beta-parameters) given the five global routing parameters (gamma) at the level 0 scale.
!> \authors Luis Samaniego, Stephan Thober
!> \date Aug 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
module mo_mrm_mpr
  use mo_kind, only : dp
  use mo_message, only : message

  implicit none

  public :: reg_rout
  public :: mrm_init_param
  public :: mrm_update_param
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
  !>       \param[in] "real(dp), dimension(:) :: fFPimp" fraction of the flood plain with
  !>       impervious layer
  !>       \param[in] "real(dp) :: TS"                   - [h] time step in

  !    INTENT(OUT)
  !>       \param[out] "real(dp), dimension(:) :: C1" routing parameter C1 (Chow, 25-41)
  !>       \param[out] "real(dp), dimension(:) :: C2" routing parameter C2 (")

  !    HISTORY
  !>       \authors Stephan Thober, Rohini Kumar

  !>       \date Dec 2012

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting

  subroutine reg_rout(param, length, slope, fFPimp, TS, C1, C2)
    implicit none

    ! input parameter
    real(dp), dimension(5), intent(in) :: param

    ! [m] total length
    real(dp), dimension(:), intent(in) :: length

    ! average slope
    real(dp), dimension(:), intent(in) :: slope

    ! fraction of the flood plain with
    ! impervious layer
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

  ! --------------------------------------------------------------------------
  ! L11 PARAMETERS
  ! --------------------------------------------------------------------------
  ! The parameters are set following Thober et al. 2017
  ! Modified:
  !    NAME
  !        mrm_init_param

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"           domain number
  !>       \param[in] "real(dp), dimension(:) :: param" input parameter (param(1) is celerity in m/s)

  !    HISTORY
  !>       \authors Stephan Thober, Matthias Kelbling

  !>       \date Jun 2018

  ! Modifications:

  subroutine mrm_init_param(iDomain, param)

    use mo_constants, only : HourSecs
    use mo_common_mHM_mRM_variables, only : resolutionRouting, timeStep, optimize
    use mo_common_variables, only : iFlag_cordinate_sys, domainMeta, processMatrix
    use mo_kind, only : dp, i4
    use mo_mrm_constants, only : given_TS
    use mo_mrm_global_variables, only : level11, L11_tsRout, domain_mrm, L11_celerity
    use mo_string_utils, only : num2str
    use mo_utils, only : locate, notequal
    use mo_mrm_net_startup, only : L11_calc_celerity

    implicit none

    ! domain number
    integer(i4), intent(in) :: iDomain

    ! input parameter (param(1) is celerity in m/s)
    real(dp), dimension(:), intent(in) :: param

    ! index selected from given_TS
    integer(i4) :: index

    ! spatial routing resolution
    real(dp) :: deltaX

    ! [s] wave travel time parameter
    real(dp) :: K

    ! start and end index at level11
    integer(i4) :: s11, e11

    ! Number of cells within

    ! initialize indices
    s11 = level11(iDomain)%iStart
    e11 = level11(iDomain)%iEnd

    ! temporal resolution of routing
    if (iDomain .eq. 1 .and. .not. allocated(L11_tsRout)) then
      allocate(L11_tsRout(domainMeta%nDomains))
      L11_TSrout = 0._dp
    end if

    if (processMatrix(8, 1) .eq. 1) then
       L11_tsRout = timestep * HourSecs

       if ( NOTEQUAL(mod(HourSecs * 24.0_dp, L11_tsRout(iDomain)), 0.0_dp) .and. &
            (domain_mrm(iDomain)%nInflowGauges .gt. 0)) then
          call message('***WARNING: routing timestep is not a multiple of 24 h.')
          call message('            Inflowgauge timeseries is averaged over values')
          call message('            of different days, small mismatches at')
          call message('            inflowgauge location might occur.')
       end if

    else

      ! called for initialization
      call mrm_update_param(iDomain, param)

    end if

    call message('')
    call message('    domain: '//num2str(iDomain, '(i3)'))
    call message('      routing resolution [s]:. '//num2str(L11_tsRout(iDomain), '(f7.0)'))
    call message('      routing factor:......... '//num2str(L11_tsRout(iDomain) / (timestep * HourSecs), '(f5.2)'))

    if ( NOTEQUAL(mod(HourSecs * 24.0_dp, L11_tsRout(iDomain)), 0.0_dp) .and. &
        (domain_mrm(iDomain)%nInflowGauges .gt. 0)) then
       call message('***WARNING: routing timestep is not a multiple of 24 h.')
       call message('            Inflowgauge timeseries is averaged over values')
       call message('            of different days, small mismatches at')
       call message('            inflowgauge location might occur.')
    end if

  end subroutine mrm_init_param

  !    NAME
  !        mrm_update_param

  !    PURPOSE
  !>       \brief TODO: add description

  !>       \details TODO: add description

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain"           domain number
  !>       \param[in] "real(dp), dimension(1) :: param" celerity parameter [m s-1]

  !    HISTORY
  !>       \authors Stehpan Thober, Matthias Kelbling

  !>       \date Jun 2018

  ! Modifications:

  subroutine mrm_update_param(iDomain, param)

    use mo_kind, only: i4, dp
    use mo_common_variables, only: processMatrix, iFlag_cordinate_sys
    use mo_mrm_global_variables, only: &
         ! input variable
         level11, &
         L11_TSrout, &
         L11_celerity, L11_nOutlets, L11_length, &
         ! output variables
         L11_C1, L11_C2
    use mo_common_mHM_mRM_variables, only: resolutionRouting, optimize, timeStep, &
         optimize
    use mo_mrm_constants, only: rout_space_weight, given_TS
    use mo_utils, only: locate
    use mo_mrm_net_startup, only: L11_calc_celerity
    use mo_mrm_constants, only: given_TS
    use mo_constants, only: HourSecs
    use mo_string_utils, only: num2str
    use mo_utils, only: locate

    implicit none

    ! domain number
    integer(i4), intent(in) :: iDomain

    ! celerity parameter [m s-1]
    real(dp), intent(in), dimension(1) :: param

    integer(i4) :: s11
    integer(i4) :: e11
    integer(i4) :: nNodes

    ! index selected from given_TS
    integer(i4) :: ind

    ! spatial routing resolution
    real(dp) :: deltaX
    real(dp), allocatable :: length(:)

    ! [s] wave travel time parameter
    real(dp), allocatable :: K(:)

    ! [1] Muskingum diffusion parameter (attenuation)
    real(dp) :: xi

    ! get domain information
    s11 = level11(iDomain)%iStart
    e11 = level11(iDomain)%iEnd
    Nnodes = level11(iDomain)%nCells

    allocate(K(nNodes))

    if (ProcessMatrix(8, 1) .eq. 2_i4) then

      ! [s] wave travel time parameter
      K(:) = L11_length(s11: e11) / param(1)

    else if (ProcessMatrix(8, 1) .eq. 3_i4) then

      ! [s] wave travel time parameter
      call L11_calc_celerity( iDomain, param)

      ! Allocate and calculate K
      K(:) = L11_length(s11: e11) / L11_celerity(s11:e11)

    end if

    ! set time-weighting scheme
    xi = abs(rout_space_weight) ! set weighting factor to 0._dp

    ! determine routing timestep
    ind = locate(given_TS, minval(K(1:(nNodes-L11_nOutlets(iDomain)))))

    ! set min-wave traveltime to min given_TS
    if (ind .lt. 1) ind = 1
    L11_TSrout(iDomain) = given_TS(ind)

    ! Muskingum parameters
    L11_C1(s11:e11) = L11_TSrout(iDomain) / ( K(:) * (1.0_dp - xi) + 0.5_dp * L11_TSrout(iDomain) )
    L11_C2(s11:e11) = 1.0_dp - L11_C1(s11:e11) * K(:) / L11_TSrout(iDomain)

    deallocate(K)

    ! optional print
    ! print *, 'C1 Muskingum routing parameter: ', L11_C1(s11)
    ! print *, 'C2 Muskingum routing parameter: ', L11_C2(s11)

  end subroutine mrm_update_param

end module mo_mrm_mpr
