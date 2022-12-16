!> \file    mo_mhm_bfi.f90
!> \brief   \copybrief mo_mhm_bfi
!> \details \copydetails mo_mhm_bfi

!> \brief   Module to calculate BFI form gauging stations in mHM.
!> \version 0.1
!> \authors Sebastian Mueller
!> \date    Apr 2022
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
module mo_mhm_bfi

  use mo_kind,       only: i4, dp

  implicit none
  private
  public :: calculate_BFI

contains

  !> \brief Calculate BFI from given discharge observation.
  subroutine calculate_BFI()
    use mo_common_mhm_mrm_variables, only : evalPer, nTstepDay, warmingDays
    use mo_common_variables, only : domainMeta
    use mo_global_variables, only : BFI_obs
    use mo_message, only : error_message, message
    use mo_string_utils, only : num2str
    use mo_utils, only : ge
    use mo_mrm_global_variables, only : gauge, nMeasPerDay
    use mo_orderpack, only: sort_index
    use mo_eckhardt_filter, only: eckhardt_filter_fit, BFI

    implicit none
    real(dp), allocatable :: baseflow(:), runoff_obs(:)
    logical, allocatable :: runoff_obs_mask(:)
    integer(i4), allocatable :: id_sort(:)
    integer(i4) :: iDomain, id, length, mask_length, tt

    ! check for daily timesteps
    if ( nMeasPerDay > 1 ) call error_message("BFI: calculation only possible with daily discharge values!")

    ! extract domain Id from gauge Id
    if (size(gauge%domainId) /= domainMeta%nDomains) call error_message("BFI: number of gauges and domains need to be equal!")
    allocate(id_sort(size(gauge%domainId)))
    id_sort = sort_index(gauge%domainId)

    call message()
    do iDomain = 1, domainMeta%nDomains
      ! skip current calculation if BFI is given
      if ( .not. (BFI_obs(iDomain) < 0.0_dp) ) then
        call message( &
          "  BFI: using given BFI value ", &
          trim(adjustl(num2str(BFI_obs(iDomain)))), &
          " for domain ", &
          trim(adjustl(num2str(iDomain))) &
        )
        cycle
      end if

      id = id_sort(iDomain)
      if (gauge%domainId(id) /= iDomain) call error_message("BFI: exactly one gauge per domain required!")

      ! get length of evaluation period times
      length = (evalPer(iDomain)%julEnd - evalPer(iDomain)%julStart + 1) * nMeasPerDay

      ! extract measurements
      if (allocated(runoff_obs_mask)) deallocate(runoff_obs_mask)
      if (allocated(runoff_obs)) deallocate(runoff_obs)
      if (allocated(baseflow)) deallocate(baseflow)
      allocate(runoff_obs_mask(length))
      allocate(runoff_obs(length))
      runoff_obs = gauge%Q(1 : length, id)

      ! create mask of observed runoff
      forall(tt = 1 : length) runoff_obs_mask(tt) = ge(runoff_obs(tt), 0.0_dp)

      ! calculate BFI
      baseflow = eckhardt_filter_fit(runoff_obs, mask=runoff_obs_mask)
      BFI_obs(iDomain) = BFI(baseflow, runoff_obs, mask=runoff_obs_mask)

      call message( &
        "  BFI:  calculated BFI value ", &
        trim(adjustl(num2str(BFI_obs(iDomain)))), &
        " for domain ", &
        trim(adjustl(num2str(iDomain))) &
      )
    end do

  end subroutine calculate_BFI

end module mo_mhm_bfi
