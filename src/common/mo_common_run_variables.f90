!> \file mo_common_run_variables.f90
!> \brief \copybrief mo_common_run_variables
!> \details \copydetails mo_common_run_variables

!> \brief Provides structures needed by mhm_eval to store current run config.
!> \author Sebastian Mueller
!> \date Jan 2022
!> \version 0.1
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_common_run_variables

  use mo_kind, only : i4, dp
  use mo_common_datetime_type, only : datetimeinfo
  use mo_nc_output, only : OutputDataset

  implicit none

  !> \class   run_cfg_t
  !> \brief   This is a container to hold all information while running mHM.
  type run_cfg_t
    !> time step counter
    integer(i4) :: time_step
    !> currently selected domain
    integer(i4) :: selected_domain
    !> number of domains simulated in this mhm_eval run. Depends on opti_function
    integer(i4) :: nDomains
    !> output runoff
    logical :: output_runoff = .false.
    !> output BFI
    logical :: output_BFI = .false.
    !> currently used parameter set
    real(dp), dimension(:), allocatable :: parameterset
    !> selected domains
    integer(i4), dimension(:), allocatable :: domain_indices
    !> fraction of NOT sealed area
    real(dp), dimension(:, :, :), allocatable :: L1_fNotSealed
    !> output mHM NetCDF object
    type(OutputDataset) :: nc_mhm
    !> output mRM NetCDF object
    type(OutputDataset) :: nc_mrm
    !> output groundwater NetCDF object
    type(OutputDataset) :: nc_gw
    !> No. of cells at level 1 for current Domain
    integer(i4) :: nCells
    !> start and end index at level 1 for current Domain
    integer(i4) :: s1, e1
    !> meteorological time step for process 5 (PET)
    integer(i4), dimension(6) :: iMeteo_p5
    !> process 5: start and end index of vectors
    !! index 1: pet
    !! index 2: tmin
    !! index 3: tmax
    !! index 4: netrad
    !! index 5: absolute vapour pressure
    !! index 6: windspeed
    integer(i4), dimension(6) :: s_p5, e_p5
    !> start and end index of meteo variables
    integer(i4) :: s_meteo, e_meteo
    !> pointer to current domain L1 mask
    logical, dimension(:, :), pointer :: mask1
    !> index of meteo time-step
    integer(i4) :: iMeteoTS
    !> datetimeinfo variable for everything that has to do with time dependend calculations
    type(datetimeinfo) :: domainDateTime
    !> discharge timestep
    integer(i4) :: iDischargeTS
    !> start and end index at L11
    integer(i4) :: s11, e11
    !> factor between routing and hydrological modelling resolution
    real(dp) :: tsRoutFactor
    !> factor between routing and hydrological modelling resolution (dummy)
    real(dp) :: tsRoutFactorIn
    !> timestep of runoff to rout [h]
    !! - identical to timestep of input if tsRoutFactor is less than 1
    !! - tsRoutFactor * timestep if tsRoutFactor is greater than 1
    integer(i4) :: timestep_rout
    !> Runoff that is input for routing
    real(dp), allocatable, dimension(:) :: RunToRout
    !> inflowing discharge
    real(dp), allocatable, dimension(:) :: InflowDischarge
    !> pointer to current domain L11 mask
    logical, pointer, dimension(:, :) :: mask11
    !> flag for performing routing
    logical :: doRoute
  contains
    !> \copydoc mo_common_run_variables::get_domain_index
    procedure :: get_domain_index!< \see mo_common_run_variables::get_domain_index
    !> \copydoc mo_common_run_variables::clean_up
    procedure :: clean_up!< \see mo_common_run_variables::clean_up
  end type run_cfg_t

  !> This is a container to hold all information while running mHM
  type(run_cfg_t), public :: run_cfg

contains

  !> \brief get domain index from domain loop counter
  !> \return index
  integer(i4) function get_domain_index(self, i) result(idx)
    implicit none
    class(run_cfg_t), intent(in) :: self
    integer(i4), intent(in) :: i !< domain loop counter
    idx = self%domain_indices(i)
  end function get_domain_index

  !> \brief clean up run variables
  subroutine clean_up(self)
    implicit none
    class(run_cfg_t), intent(inout) :: self
    if ( allocated(self%parameterset) ) deallocate(self%parameterset)
    if ( allocated(self%domain_indices) ) deallocate(self%domain_indices)
    if ( allocated(self%L1_fNotSealed) ) deallocate(self%L1_fNotSealed)
    if ( allocated(self%RunToRout) ) deallocate(self%RunToRout)
    if ( allocated(self%InflowDischarge) ) deallocate(self%InflowDischarge)
  end subroutine clean_up

end module mo_common_run_variables
