!> \file    mo_common_mpi_tools.f90
!> \brief   \copybrief mo_common_mpi_tools
!> \details \copydetails mo_common_mpi_tools

!> \brief   tools for MPI communication that are mHM or mRM specific
!> \author  Maren Kaluza
!> \author  Sebastian Mueller
!> \date    2019-2021
!> \details This module contains sending and receiving subroutines for
!!          data that are specific for mHM or mRM
!> \copyright Copyright 2005-\today, the mHM Developers, Luis Samaniego, Sabine Attinger: All rights reserved.
!! mHM is released under the LGPLv3+ license \license_note
!> \ingroup f_common
MODULE mo_common_mpi_tools

#ifdef MPI
  use mo_kind, only : i4, dp
  use mo_message, only : message
  use mo_string_utils, only : num2str
  use mo_common_variables, only: comm
  use mpi_f08
#endif

  IMPLICIT NONE

  PRIVATE

#ifdef MPI
  PUBLIC :: distribute_parameterset, get_parameterset
#endif
  public :: mpi_tools_init
  public :: mpi_tools_finalize

  ! ------------------------------------------------------------------

contains
#ifdef MPI
  !> \brief Distrubute parameter set with MPI.
  subroutine distribute_parameterset(parameterset)
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

  !> \brief Get distrubuted parameter set with MPI.
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

  !> \brief Finalize the MPI run of mHM.
  subroutine mpi_tools_init()

    implicit none

#ifdef MPI
    integer             :: ierror
    integer(i4)         :: nproc, rank

    ! Initialize MPI
    call MPI_Init(ierror)
    call MPI_Comm_dup(MPI_COMM_WORLD, comm, ierror)
    ! find number of processes nproc
    call MPI_Comm_size(comm, nproc, ierror)
    ! find the number the process is referred to, called rank
    call MPI_Comm_rank(comm, rank, ierror)
    call message('MPI!, comm ', num2str(rank), num2str(nproc))
#endif

  end subroutine mpi_tools_init

  !> \brief Finalize the MPI run of mHM.
  subroutine mpi_tools_finalize()

    implicit none

#ifdef MPI
    integer             :: ierror
    integer(i4)         :: nproc, rank

    ! find number of processes nproc
    call MPI_Comm_size(comm, nproc, ierror)
    call MPI_Comm_rank(comm, rank, ierror)
    call message('MPI finished ', num2str(rank), num2str(nproc))
    call MPI_Finalize(ierror)
#endif

  end subroutine mpi_tools_finalize

END MODULE mo_common_mpi_tools
