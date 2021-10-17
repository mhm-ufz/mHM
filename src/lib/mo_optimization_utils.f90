module mo_optimization_utils

  use mo_kind, only : dp

  implicit none

  abstract interface
    subroutine eval_interface(parameterset, opti_domain_indices, runoff, smOptiSim, neutronsOptiSim, etOptiSim, twsOptiSim)
      use mo_kind, only : dp, i4
      use mo_optimization_types, only : optidata_sim
      real(dp),    dimension(:), intent(in) :: parameterset
      integer(i4), dimension(:),                 optional, intent(in)  :: opti_domain_indices
      real(dp),    dimension(:, :), allocatable, optional, intent(out) :: runoff        ! dim1=time dim2=gauge
      type(optidata_sim), dimension(:), optional, intent(inout) :: smOptiSim       ! dim1=ncells, dim2=time
      type(optidata_sim), dimension(:), optional, intent(inout) :: neutronsOptiSim ! dim1=ncells, dim2=time
      type(optidata_sim), dimension(:), optional, intent(inout) :: etOptiSim       ! dim1=ncells, dim2=time
      type(optidata_sim), dimension(:), optional, intent(inout) :: twsOptiSim      ! dim1=ncells, dim2=time
    end subroutine
  end interface

  interface
    function objective_interface (parameterset, eval, arg1, arg2, arg3)
      use mo_kind, only : dp
      import eval_interface
      real(dp), intent(in), dimension(:) :: parameterset
      procedure(eval_interface), INTENT(IN), pointer :: eval
      real(dp), optional, intent(in) :: arg1
      real(dp), optional, intent(out) :: arg2
      real(dp), optional, intent(out) :: arg3

      real(dp) :: objective_interface
    end function objective_interface
  end interface

end module mo_optimization_utils


