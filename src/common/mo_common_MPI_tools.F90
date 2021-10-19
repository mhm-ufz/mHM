!>       \file use mo_common_mpi_tools.f90

!>       \brief tools for MPI communication that are mHM or mRM specific

!>       \details This module contains sending and receiving subroutines for
!>                data that are specific for mHM or mRM

!>       \authors Maren Kaluza

!>       \date Jul 2019

! Modifications:

MODULE mo_common_mpi_tools

  use mo_kind, only : i4, dp

  IMPLICIT NONE

  PRIVATE

#ifdef MPI
  PUBLIC :: distribute_parameterset, get_parameterset
#endif

  ! ------------------------------------------------------------------

contains
#ifdef MPI
  subroutine distribute_parameterset(parameterset)
    use mpi_f08
    use mo_common_variables, only : domainMeta
    real(dp), dimension(:),    intent(in) :: parameterset

    integer(i4) :: nproc, iproc, dimen
    integer(i4) :: ierror

    call MPI_Comm_size(domainMeta%comMaster, nproc, ierror)
    dimen = size(parameterset(:))
    do iproc = 1, nproc-1
      call MPI_Send(dimen, 1, &
                    MPI_INTEGER,iproc,0,domainMeta%comMaster,ierror)
      call MPI_Send(parameterset(:),dimen, &
                    MPI_DOUBLE_PRECISION,iproc,0,domainMeta%comMaster,ierror)
    end do
  end subroutine distribute_parameterset

  subroutine get_parameterset(parameterset)
    use mpi_f08
    use mo_common_variables, only : domainMeta
    real(dp), dimension(:), allocatable, intent(inout) :: parameterset

    integer(i4) :: dimen
    integer(i4) :: ierror
    type(MPI_Status) :: status

    call MPI_Recv(dimen, 1, MPI_INTEGER, 0, 0, domainMeta%comMaster, status, ierror)
    allocate(parameterset(dimen))
    call MPI_Recv(parameterset, dimen, MPI_DOUBLE_PRECISION, 0, 0, domainMeta%comMaster, status, ierror)
  end subroutine get_parameterset
#endif
END MODULE mo_common_mpi_tools
