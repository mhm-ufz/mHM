!> \file    mo_mrm_river_head.f90
!> \brief   \copybrief mo_mrm_river_head
!> \details \copydetails mo_mrm_river_head

!> \brief   River head calculation
!> \details Enables river - groundwater interaction in mRM.
!> \authors Lennart Schueler
!> \date    Jul 2018
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mrm
module mo_mrm_river_head
  use mo_common_variables,     only : level0, domainMeta
  use mo_mrm_global_variables, only : L0_L11_remap, L11_bankfull_runoff_in, &
                                      L0_slope, L0_channel_elevation
  use mo_kind,                 only : i4, dp
  use mo_common_constants,     only : nodata_dp
  use mo_append,               only : append

  implicit none

  private

  public ::  init_masked_zeros_l0
  public :: calc_channel_elevation
  public :: calc_river_head

  contains

  !> \brief allocates memory for L0 variable
  subroutine init_masked_zeros_l0(iDomain, data)
    integer(i4), intent(in) :: iDomain
    real(dp), dimension(:), allocatable, intent(inout) :: data
    real(dp), dimension(:), allocatable :: dummy_1D

    allocate(dummy_1D(level0(iDomain)%nCells))
    dummy_1D = nodata_dp
    call append(data, dummy_1D)
    deallocate(dummy_1D)
  end subroutine init_masked_zeros_l0

  !> \brief calculates the channel elevation from the bankfull river discharge
  subroutine calc_channel_elevation()
    use mo_common_constants, only : nodata_i4
    use mo_common_variables, only : domainMeta, L0_elev
    use mo_mrm_global_variables, only : L0_fDir, L0_fAcc, L0_channel_depth
    real(dp), dimension(:,:), allocatable :: channel_dpth
    real(dp), dimension(:,:), allocatable :: channel_elev
    real(dp), dimension(:,:), allocatable :: slope
    real(dp), dimension(:,:), allocatable :: elev0
    integer(i4), dimension(:,:), allocatable :: fDir0
    integer(i4), dimension(:,:), allocatable :: fAcc0
    real(dp) n ! Manning's roughness coefficient
    integer(i4) :: nrows0, ncols0
    integer(i4) :: s0, e0
    integer(i4) i, j, k
    integer(i4) iDomain

    n = .045_dp ! m^-1/3 s from Sutanudjaja et al. 2011

    do iDomain = 1, domainMeta%nDomains
      nrows0 = level0(iDomain)%nrows
      ncols0 = level0(iDomain)%ncols
      s0 = level0(iDomain)%iStart
      e0 = level0(iDomain)%iEnd
      allocate(channel_dpth(nrows0, ncols0))
      allocate(channel_elev(nrows0, ncols0))
      allocate(elev0(nrows0, ncols0))
      allocate(fDir0(nrows0, ncols0))
      allocate(fAcc0(nrows0, ncols0))
      allocate(slope(nrows0, ncols0))
      channel_dpth(:,:) = nodata_dp
      channel_elev(:,:) = nodata_dp
      slope(:,:) = nodata_dp

      elev0(:,:) = unpack(L0_elev(s0:e0), level0(iDomain)%mask, nodata_dp)
      fDir0(:,:) = unpack(L0_fDir(s0:e0), level0(iDomain)%mask, nodata_i4)
      fAcc0(:,:) = unpack(L0_fAcc(s0:e0), level0(iDomain)%mask, nodata_i4)

      do k = 1, level0(iDomain)%nCells
        i = level0(iDomain)%CellCoor(k, 1)
        j = level0(iDomain)%CellCoor(k, 2)
        if (fAcc0(i,j) > 1) then
          slope(i,j) = calc_slope(iDomain, elev0, fDir0, i, j)
          channel_dpth(i,j) = ((n * &
  L11_bankfull_runoff_in(L0_L11_remap(iDomain)%lowres_id_on_highres(i,j))) / &
  (4.8 * slope(i,j)**.5))**.6
          channel_elev(i,j) = elev0(i,j) - channel_dpth(i,j)
        end if
      end do

      call append(L0_channel_depth, pack(channel_dpth(:,:), &
                  level0(iDomain)%mask))
      call append(L0_channel_elevation, pack(channel_elev(:,:), &
                  level0(iDomain)%mask))
      call append(L0_slope, pack(slope(:,:), &
                  level0(iDomain)%mask))

      deallocate(channel_dpth)
      deallocate(channel_elev)
      deallocate(elev0)
      deallocate(fDir0)
      deallocate(fAcc0)
      deallocate(slope)
    end do
  end subroutine calc_channel_elevation

  !> \brief calculates the river head
  subroutine calc_river_head(iDomain, L11_Qmod, river_head)
    integer(i4), intent(in) :: iDomain
    real(dp), dimension(:), intent(in) :: L11_Qmod
    real(dp), dimension(:), allocatable, intent(inout) :: river_head
    real(dp) :: n ! Manning's roughness coefficient
    integer(i4) :: s0, e0
    integer(i4) i, j, k, L11_ind

    n = .045_dp ! m^-1/3 s from Sutanudjaja et al. 2011

    s0 = level0(iDomain)%iStart
    e0 = level0(iDomain)%iEnd

    do k = s0, e0
      i = level0(iDomain)%CellCoor(k - s0 + 1, 1)
      j = level0(iDomain)%CellCoor(k - s0 + 1, 2)
      if (i >= 0 .and. i < 99999 .and. j  >= 0 .and. j < 99999) then
        L11_ind = L0_L11_remap(iDomain)%lowres_id_on_highres(i,j)
        ! TODO L11_Qmid(L11_ind) causes IEEE_UNDERFLOW_FLAG IEEE_DENORMAL
        river_head(k) = L0_channel_elevation(k) + &
            (n * L11_Qmod(L11_ind) / L11_bankfull_runoff_in(L11_ind) / &
            L0_slope(k)**.5)**.6
      end if
    end do

  end subroutine calc_river_head

  !> \brief calculates domain slope
  function calc_slope(iDomain, elev0, fDir0, i, j) result(slope)
    use mo_common_variables, only: iFlag_cordinate_sys
    use mo_mrm_net_startup,      only: cellLength, moveDownOneCell
    integer(i4), intent(in) :: iDomain
    integer(i4), intent(in) :: i, j
    integer(i4), intent(in), dimension(:,:), allocatable :: fDir0
    real(dp), intent(in), dimension(:,:), allocatable :: elev0
    real(dp) :: slope, length
    integer(i4) :: i_down, j_down

    call cellLength(iDomain, fDir0(i, j), i, j, &
                    iFlag_cordinate_sys, length)
    i_down = i
    j_down = j
    call moveDownOneCell(fDir0(i, j), i_down, j_down)

    slope = (elev0(i,j) - elev0(i_down, j_down)) / length
    ! TODO: as soon as current gfortran compiler is available on EVE,
    ! use ieee_isnan from ieee_arithmetic module, instead of
    ! slope /= slope
    if(slope < 0.0001_dp .OR. slope > 20._dp .OR. slope /= slope) then
      slope = 0.0001_dp
    end if
  end function calc_slope

end module mo_mrm_river_head
