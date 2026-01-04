!> \file mo_optimization_utils.f90
!> \copydoc mo_optimization_utils

!> \brief Utility functions, such as interface definitions, for optimization routines.
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
module mo_optimization_utils

  use mo_kind, only : dp
  use mo_optimizee, only: optimizee
  use mo_message, only : error_message

  implicit none

  !> \brief Interface for evaluation routine.
  abstract interface
    subroutine eval_interface(parameterset, opti_domain_indices, runoff, smOptiSim, neutronsOptiSim, etOptiSim, twsOptiSim, BFI, &
      sweOptiSim)
      use mo_kind, only : dp, i4
      use mo_optimization_types, only : optidata_sim
      real(dp),    dimension(:), intent(in) :: parameterset
      integer(i4), dimension(:),                 optional, intent(in)  :: opti_domain_indices
      real(dp),    dimension(:, :), allocatable, optional, intent(out) :: runoff   !< dim1=time dim2=gauge
      type(optidata_sim), dimension(:), optional, intent(inout) :: smOptiSim       !< dim1=ncells, dim2=time
      type(optidata_sim), dimension(:), optional, intent(inout) :: neutronsOptiSim !< dim1=ncells, dim2=time
      type(optidata_sim), dimension(:), optional, intent(inout) :: etOptiSim       !< dim1=ncells, dim2=time
      type(optidata_sim), dimension(:), optional, intent(inout) :: twsOptiSim      !< dim1=ncells, dim2=time
      type(optidata_sim), dimension(:), optional, intent(inout) :: sweOptiSim      !< dim1=ncells, dim2=time
      real(dp),    dimension(:), allocatable, optional, intent(out) :: BFI         !< baseflow index, dim1=domainID
    end subroutine
  end interface

  !> \brief Interface for objective function.
  interface
    function objective_interface (parameterset, eval, arg1, arg2, arg3)
      use mo_kind, only : dp
      import eval_interface
      real(dp), intent(in), dimension(:) :: parameterset !< parameter set
      procedure(eval_interface), INTENT(IN), pointer :: eval !< evaluation routine
      real(dp), optional, intent(in) :: arg1 !< optional argument 1
      real(dp), optional, intent(out) :: arg2 !< optional argument 2
      real(dp), optional, intent(out) :: arg3 !< optional argument 3

      real(dp) :: objective_interface
    end function objective_interface
  end interface

  !> \brief Optimizee for a eval-objective pair
  type, extends(optimizee) :: mhm_optimizee
    procedure(eval_interface), pointer, nopass :: eval_pointer => null()  !< Pointer to the eval
    procedure(objective_interface), pointer, nopass :: obj_pointer => null()  !< Pointer to the objective
  contains
    procedure :: evaluate => evaluate_obj_eval
  end type mhm_optimizee

  contains

  !> \brief Implementation of the evaluate procedure for a eval-objective pair
  function evaluate_obj_eval(self, parameters, sigma, stddev_new, likeli_new) result(value)
    class(mhm_optimizee), intent(inout) :: self
    real(DP), dimension(:), intent(in) :: parameters
    real(DP), intent(in),  optional :: sigma
    real(DP), intent(out), optional :: stddev_new
    real(DP), intent(out), optional :: likeli_new
    real(DP) :: value

    ! Ensure the eval function pointer is set
    if (.not. associated(self%eval_pointer)) then
      call error_message("Eval function pointer is not set in mhm_optimizee!")
    end if
    ! Ensure the objective function pointer is set
    if (.not. associated(self%obj_pointer)) then
      call error_message("Objective function pointer is not set in mhm_optimizee!")
    end if

    ! Call the objective function pointer
    value = self%obj_pointer(parameters, self%eval_pointer, sigma, stddev_new, likeli_new)
  end function evaluate_obj_eval

end module mo_optimization_utils
