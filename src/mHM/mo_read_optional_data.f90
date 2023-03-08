!> \file mo_read_optional_data.f90
!> \brief   \copybrief mo_read_optional_data
!> \details \copydetails mo_read_optional_data

!> \brief Read optional data for mHM calibration.
!> \details Data have to be provided in resolution of the hydrology.
!> \authors Matthias Zink
!> \date Mar 2015
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_mhm
MODULE mo_read_optional_data

  USE mo_kind, ONLY : i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: readOptidataObs

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !    NAME
  !        readOptidataObs

  !    PURPOSE
  !>       \brief Read evapotranspiration data from NetCDF file for calibration

  !>       \details This routine reads oberved evapotranspiration fields which are used for model
  !>       calibration. The evapotranspiration file is expected to be called "et.nc" with
  !>       a variable "et" inside. The data are read only for the evaluation period
  !>       they are intended to be used for calibration. Evapotranspiration data are only
  !>       read if one of the corresponding objective functions is chosen.

  !    INTENT(IN)
  !>       \param[in] "integer(i4) :: iDomain" domain Id

  !    HISTORY
  !>       \authors Johannes Brenner

  !>       \date Feb 2017

  ! Modifications:
  ! Robert Schweppe Jun 2018 - refactoring and reformatting
  ! Maren Kaluza Oct 2019 - copied from evapotranspiration and adopted for tws
  ! Maren Kaluza Nov 2019 - removed redundant code, reading works for any gridded optidata

  subroutine readOptidataObs(iDomain, domainID, L1_optiObs)

    use mo_common_mhm_mrm_variables, only : evalPer
    use mo_common_variables, only : level1
    use mo_optimization_types, only: optidata
    use mo_message, only : message
    use mo_read_nc, only : read_nc
    use mo_string_utils, only : num2str
    use mo_timer, only : timer_get, timer_start, &
                         timer_stop

    implicit none

    ! domain Id
    integer(i4), intent(in) :: iDomain

    ! domain Id
    integer(i4), intent(in) :: domainID

    ! opti data
    type(optidata), intent(inout) :: L1_optiObs

    ! loop  vars packing L1_data to L1_data_packed
    integer(i4) :: t

    ! level 1 number of culomns and rows
    integer(i4) :: nrows1, ncols1

    ! mask of level 1 for packing
    logical, dimension(:, :), allocatable :: mask1

    ! ncells1 of level 1
    integer(i4) :: ncells1

    ! data at level-1
    real(dp), dimension(:, :, :), allocatable :: L1_data

    ! mask at level-1
    logical, dimension(:, :, :), allocatable :: L1_mask

    integer(i4) :: nTimeSteps_L1_opti      ! [-] number of time steps in L1_opti_mask

    ! get basic domain information at level-1
    nrows1 = level1(iDomain)%nrows
    ncols1 = level1(iDomain)%ncols
    ncells1 = level1(iDomain)%ncells
    mask1 = level1(iDomain)%mask

    !  domain characteristics and read meteo header
    call message('  Reading', trim(L1_optiObs%varname) ,'for domain:           ', trim(adjustl(num2str(domainID))), ' ...')
    call timer_start(1)
    call read_nc(L1_optiObs%dir, nRows1, nCols1, trim(L1_optiObs%varname), mask1, L1_data, &
            target_period = evalPer(iDomain), nctimestep = L1_optiObs%timeStepInput, nocheck = .TRUE., maskout = L1_mask)

    ! pack variables
    nTimeSteps_L1_opti = size(L1_data, 3)
    allocate(L1_optiObs%dataObs(nCells1, nTimeSteps_L1_opti))
    allocate(L1_optiObs%maskObs(nCells1, nTimeSteps_L1_opti))
    do t = 1, nTimeSteps_L1_opti
      L1_optiObs%dataObs(:, t) = pack(L1_data(:, :, t), MASK = mask1(:, :))
      L1_optiObs%maskObs(:, t) = pack(L1_mask(:, :, t), MASK = mask1(:, :))
    end do

    !free space
    deallocate(L1_data)

    call timer_stop(1)
    call message('    in ', trim(num2str(timer_get(1), '(F9.3)')), ' seconds.')

  end subroutine readOptidataObs

END MODULE mo_read_optional_data
