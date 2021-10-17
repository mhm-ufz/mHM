!>       \file mo_common_mhm_mrm_restart.f90

!>       \brief TODO: add description

!>       \details TODO: add description

!>       \authors Robert Schweppe

!>       \date Aug 2019

! Modifications:

module mo_common_mHM_mRM_restart

  use mo_kind, only : i4, dp
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: check_dimension_consistency
  INTERFACE check_consistency_element
    MODULE PROCEDURE check_consistency_element_i4, &
            check_consistency_element_dp
  end interface check_consistency_element


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        check_dimension_consistency

  !    PURPOSE
  !>       \brief checks dimension configurations read from restart file

  !    HISTORY
  !>       \authors Robert Schweppe

  !>       \date Aug 2019

  ! Modifications:

  subroutine check_dimension_consistency(iBasin, nSoilHorizons_temp, soilHorizonBoundaries_temp, &
          nLAIs_temp, LAIBoundaries_temp, nLandCoverPeriods_temp, landCoverPeriodBoundaries_temp)
    use mo_mpr_global_variables, only: nSoilHorizons_mHM, HorizonDepth_mHM, nLAI, LAIBoundaries
    use mo_common_variables, only: nLCoverScene, LC_year_start, LC_year_end
    use mo_string_utils, only: compress, num2str

    integer(i4), intent(in) :: iBasin

    integer(i4), intent(in) :: nSoilHorizons_temp, nLAIs_temp, nLandCoverPeriods_temp
    real(dp), dimension(:), intent(inout) :: landCoverPeriodBoundaries_temp, soilHorizonBoundaries_temp, &
            LAIBoundaries_temp
    character(256) :: errorString

    integer(i4) :: k

    ! compare local to global
    call check_consistency_element(nLCoverScene, nLandCoverPeriods_temp, 'number of land cover periods', iBasin)
    call check_consistency_element(nSoilHorizons_mHM, nSoilHorizons_temp, 'number of soil horizons', iBasin)
    call check_consistency_element(nLAI, nLAIs_temp, 'number of LAI timesteps', iBasin)

    ! now check the boundaries
    do k=1, nLCoverScene
      errorString = compress(trim(num2str(k)))//'th land cover boundary'
      call check_consistency_element(real(LC_year_start(k), dp), landCoverPeriodBoundaries_temp(k), errorString, iBasin)
    end do
    errorString = 'last land cover boundary (with 1 year added due to real/int conversion) '
    call check_consistency_element(real(LC_year_end(nLCoverScene) + 1_i4, dp), &
            landCoverPeriodBoundaries_temp(nLCoverScene+1), errorString, iBasin)

    ! last soil horizon is spatially variable, so this is not checked yet
    ! first soil horizon 0 and not contained in HorizonDepth_mHM, so skip that, too
    do k=2, nSoilHorizons_mHM
      errorString = compress(trim(num2str(k)))//'th soil horizon boundary'
      call check_consistency_element(HorizonDepth_mHM(k-1), soilHorizonBoundaries_temp(k), errorString, iBasin)
    end do

    do k=1, nLAI+1
      errorString = compress(trim(num2str(k)))//'th LAI period boundary'
      call check_consistency_element(LAIBoundaries(k), LAIBoundaries_temp(k), errorString, iBasin)
    end do


  end subroutine check_dimension_consistency

  subroutine check_consistency_element_dp(item1, item2, name, iBasin)
    use mo_utils, only: ne
    use mo_string_utils, only: compress, num2str
    use mo_message, only: message

    real(dp), intent(in) :: item1, item2
    character(*), intent(in) :: name
    integer(i4), intent(in) :: iBasin

    if (ne(item1, item2)) then
      call message('The ', trim(name),&
                  ' as set in the configuration file (', &
                  compress(trim(num2str(item1))), &
                  ') does not conform with basin ', &
                  compress(trim(num2str(iBasin))), ' (', compress(trim(num2str(item2))), ').')
      stop 1
    end if
  end subroutine check_consistency_element_dp

  subroutine check_consistency_element_i4(item1, item2, name, iBasin)
    use mo_utils, only: ne
    use mo_string_utils, only: compress, num2str
    use mo_message, only: message

    integer(i4), intent(in) :: item1, item2
    character(*), intent(in) :: name
    integer(i4), intent(in) :: iBasin

    if (item1 /= item2) then
      call message('The ', trim(name),&
                  ' as set in the configuration file (', &
                  compress(trim(num2str(item1))), &
                  ') does not conform with basin ', &
                  compress(trim(num2str(iBasin))), ' (', compress(trim(num2str(item2))), ').')
      stop 1
    end if
  end subroutine check_consistency_element_i4

end module mo_common_mHM_mRM_restart


