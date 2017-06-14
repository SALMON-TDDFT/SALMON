!
!  Copyright 2016 ARTED developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
module salmon_communication
  use mpi, only: MPI_PROC_NULL
  implicit none

  integer, public, parameter :: ROOT_PROCID = 0
  integer, public, parameter :: COMM_PROC_NULL = MPI_PROC_NULL

  ! call once
  public :: comm_init
  public :: comm_finalize

  ! p2p communication
  public :: comm_send
  public :: comm_recv
  public :: comm_exchange

  ! p2p immediate communication
  public :: comm_isend
  public :: comm_irecv
  public :: comm_wait
  public :: comm_wait_all

  ! p2p persistent communication
  !public :: comm_send_init
  !public :: comm_recv_init
  !public :: comm_start_all

  ! collective communication
  public :: comm_sync_all
  public :: comm_summation
  public :: comm_bcast
  public :: comm_allgatherv
  public :: comm_get_min
  public :: comm_get_max

  ! group (communicator)
  public :: comm_get_globalinfo
  public :: comm_get_groupinfo
  public :: comm_create_group

  ! utils
  public :: comm_is_root


  type, public :: comm_maxloc_type
    real(8) :: val
    integer :: rank
  end type

  interface comm_send
    ! 5-D array
    module procedure comm_send_array5d_real8
    module procedure comm_send_array5d_complex8
  end interface

  interface comm_recv
    ! 5-D array
    module procedure comm_recv_array5d_real8
    module procedure comm_recv_array5d_complex8
  end interface

  interface comm_exchange
    ! 3-D array
    module procedure comm_exchange_array3d_real8
    module procedure comm_exchange_array3d_complex8

    ! 5-D array
    module procedure comm_exchange_array5d_real8
    module procedure comm_exchange_array5d_complex8
  end interface

  interface comm_isend
    ! 5-D array
    module procedure comm_isend_array5d_real8
    module procedure comm_isend_array5d_complex8
  end interface

  interface comm_irecv
    ! 5-D array
    module procedure comm_irecv_array5d_real8
    module procedure comm_irecv_array5d_complex8
  end interface

  interface comm_summation
    ! scalar
    module procedure comm_summation_real8
    module procedure comm_summation_complex8

    ! 1-D array
    module procedure comm_summation_array1d_integer
    module procedure comm_summation_array1d_real8
    module procedure comm_summation_array1d_complex8

    ! 2-D array
    module procedure comm_summation_array2d_integer
    module procedure comm_summation_array2d_real8
    module procedure comm_summation_array2d_complex8

    ! 3-D array
    module procedure comm_summation_array3d_integer
    module procedure comm_summation_array3d_real8
    module procedure comm_summation_array3d_complex8

    ! 4-D array
    module procedure comm_summation_array4d_real8
    module procedure comm_summation_array4d_complex8
  end interface

  interface comm_bcast
    ! scalar
    module procedure comm_bcast_integer
    module procedure comm_bcast_real8
    module procedure comm_bcast_character
    module procedure comm_bcast_logical

    ! 1-D array
    module procedure comm_bcast_array1d_integer
    module procedure comm_bcast_array1d_real8
    module procedure comm_bcast_array1d_character

    ! 2-D array
    module procedure comm_bcast_array2d_integer
    module procedure comm_bcast_array2d_real8

    ! 3-D array
    module procedure comm_bcast_array3d_real8
    module procedure comm_bcast_array3d_complex8
  end interface

  interface comm_allgatherv
    ! 1-D array
    module procedure comm_allgatherv_array1d_real8
  end interface

  interface comm_get_min
    ! 1-D array
    module procedure comm_get_min_array1d_real8
  end interface

  interface comm_get_max
    ! scalar
    module procedure comm_get_maxloc

    ! 1-D array
    module procedure comm_get_max_array1d_real8
  end interface

  interface comm_logical_and
    ! scalar
    module procedure comm_logical_and_scalar
  end interface

  private :: get_rank, error_check

#define MPI_ERROR_CHECK(x) x; call error_check(ierr)

contains
  subroutine comm_init
    use mpi
    implicit none
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Init(ierr))
  end subroutine

  subroutine comm_finalize
    use mpi
    implicit none
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Finalize(ierr))
  end subroutine

  subroutine comm_get_globalinfo(ngid, npid, nprocs)
    use mpi
    implicit none
    integer, intent(out) :: ngid, npid, nprocs
    ngid = MPI_COMM_WORLD
    call get_rank(ngid, npid, nprocs)
  end subroutine

  subroutine comm_get_groupinfo(ngid, npid, nprocs)
    use mpi
    implicit none
    integer, intent(in)  :: ngid
    integer, intent(out) :: npid, nprocs
    call get_rank(ngid, npid, nprocs)
  end subroutine

  function comm_create_group(ngid, nprocs, key) result(ngid_dst)
    use mpi
    implicit none
    integer, intent(in) :: ngid, nprocs, key
    integer :: ngid_dst, ierr
    MPI_ERROR_CHECK(call MPI_Comm_split(ngid, nprocs, key, ngid_dst, ierr))
  end function

  function comm_is_root(npid)
    implicit none
    integer, intent(in) :: npid
    logical :: comm_is_root
    comm_is_root = npid == ROOT_PROCID
  end function

  subroutine comm_sync_all(ngid)
    use mpi
    implicit none
    integer, intent(in), optional :: ngid
    integer :: ierr
    if (present(ngid)) then
      MPI_ERROR_CHECK(call MPI_Barrier(ngid, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Barrier(MPI_COMM_WORLD, ierr))
    end if
  end subroutine


  subroutine comm_send_array5d_real8(invalue, ndest, ntag, ngroup)
    use mpi, only: MPI_REAL8, MPI_STATUS_SIZE
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Send(invalue, size(invalue), MPI_REAL8, ndest, ntag, ngroup, ierr))
  end subroutine

  subroutine comm_send_array5d_complex8(invalue, ndest, ntag, ngroup)
    use mpi, only: MPI_COMPLEX8, MPI_STATUS_SIZE
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Send(invalue, size(invalue), MPI_COMPLEX8, ndest, ntag, ngroup, ierr))
  end subroutine

  subroutine comm_recv_array5d_real8(outvalue, nsrc, ntag, ngroup)
    use mpi, only: MPI_REAL8, MPI_STATUS_SIZE
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    MPI_ERROR_CHECK(call MPI_Recv(outvalue, size(outvalue), MPI_REAL8, nsrc, ntag, ngroup, istatus, ierr))
  end subroutine

  subroutine comm_recv_array5d_complex8(outvalue, nsrc, ntag, ngroup)
    use mpi, only: MPI_COMPLEX8, MPI_STATUS_SIZE
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    MPI_ERROR_CHECK(call MPI_Recv(outvalue, size(outvalue), MPI_COMPLEX8, nsrc, ntag, ngroup, istatus, ierr))
  end subroutine

  subroutine comm_exchange_array3d_real8(invalue, nsrc, outvalue, ndest, ntag, ngroup)
    use mpi, only: MPI_REAL8, MPI_STATUS_SIZE
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ndest, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    call MPI_Sendrecv(invalue,  size(invalue),  MPI_REAL8, ndest, ntag, &
                      outvalue, size(outvalue), MPI_REAL8, nsrc,  ntag, &
                      ngroup, istatus, ierr)
    call error_check(ierr)
  end subroutine

  subroutine comm_exchange_array3d_complex8(invalue, nsrc, outvalue, ndest, ntag, ngroup)
    use mpi, only: MPI_COMPLEX8, MPI_STATUS_SIZE
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ndest, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    call MPI_Sendrecv(invalue,  size(invalue),  MPI_COMPLEX8, ndest, ntag, &
                      outvalue, size(outvalue), MPI_COMPLEX8, nsrc,  ntag, &
                      ngroup, istatus, ierr)
    call error_check(ierr)
  end subroutine

  subroutine comm_exchange_array5d_real8(invalue, nsrc, outvalue, ndest, ntag, ngroup)
    use mpi, only: MPI_REAL8, MPI_STATUS_SIZE
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ndest, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    call MPI_Sendrecv(invalue,  size(invalue),  MPI_REAL8, ndest, ntag, &
                      outvalue, size(outvalue), MPI_REAL8, nsrc,  ntag, &
                      ngroup, istatus, ierr)
    call error_check(ierr)
  end subroutine

  subroutine comm_exchange_array5d_complex8(invalue, nsrc, outvalue, ndest, ntag, ngroup)
    use mpi, only: MPI_COMPLEX8, MPI_STATUS_SIZE
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ndest, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    call MPI_Sendrecv(invalue,  size(invalue),  MPI_COMPLEX8, ndest, ntag, &
                      outvalue, size(outvalue), MPI_COMPLEX8, nsrc,  ntag, &
                      ngroup, istatus, ierr)
    call error_check(ierr)
  end subroutine

  function comm_isend_array5d_real8(invalue, ndest, ntag, ngroup) result(req)
    use mpi, only: MPI_REAL8
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: ierr, istatus, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Isend(invalue, size(invalue), MPI_REAL8, ndest, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, istatus, ierr))
  end function

  function comm_isend_array5d_complex8(invalue, ndest, ntag, ngroup) result(req)
    use mpi, only: MPI_COMPLEX8
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: ierr, istatus, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Isend(invalue, size(invalue), MPI_COMPLEX8, ndest, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, istatus, ierr))
  end function

  function comm_irecv_array5d_real8(outvalue, nsrc, ntag, ngroup) result(req)
    use mpi, only: MPI_REAL8
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: ierr, istatus, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Irecv(outvalue, size(outvalue), MPI_REAL8, nsrc, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, istatus, ierr))
  end function

  function comm_irecv_array5d_complex8(outvalue, nsrc, ntag, ngroup) result(req)
    use mpi, only: MPI_COMPLEX8
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: ierr, istatus, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Irecv(outvalue, size(outvalue), MPI_COMPLEX8, nsrc, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, istatus, ierr))
  end function

  subroutine comm_wait(req)
    use mpi, only: MPI_STATUSES_IGNORE
    implicit none
    integer, intent(in) :: req
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Wait(req, MPI_STATUSES_IGNORE, ierr))
  end subroutine

  subroutine comm_wait_all(reqs)
    use mpi, only: MPI_STATUSES_IGNORE
    implicit none
    integer, intent(in) :: reqs(:)
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Waitall(size(reqs), reqs, MPI_STATUSES_IGNORE, ierr))
  end subroutine


  subroutine comm_summation_real8(invalue, outvalue, ngroup, dest)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue
    real(8), intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, 1, MPI_REAL8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_REAL8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_complex8(invalue, outvalue, ngroup, dest)
    use mpi
    implicit none
    complex(8), intent(in)  :: invalue
    complex(8), intent(out) :: outvalue
    integer, intent(in)     :: ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, 1, MPI_COMPLEX8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_COMPLEX8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array1d_integer(invalue, outvalue, N, ngroup, dest)
    use mpi, only: MPI_INTEGER, MPI_SUM
    implicit none
    integer, intent(in)  :: invalue(:)
    integer, intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_INTEGER, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_INTEGER, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array1d_real8(invalue, outvalue, N, ngroup, dest)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array1d_complex8(invalue, outvalue, N, ngroup, dest)
    use mpi, only: MPI_COMPLEX8, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:)
    complex(8), intent(out) :: outvalue(:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_COMPLEX8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_COMPLEX8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array2d_integer(invalue, outvalue, N, ngroup, dest)
    use mpi, only: MPI_INTEGER, MPI_SUM
    implicit none
    integer, intent(in)  :: invalue(:,:)
    integer, intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_INTEGER, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_INTEGER, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array2d_real8(invalue, outvalue, N, ngroup, dest)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:,:)
    real(8), intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array2d_complex8(invalue, outvalue, N, ngroup, dest)
    use mpi, only: MPI_COMPLEX8, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:,:)
    complex(8), intent(out) :: outvalue(:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_COMPLEX8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_COMPLEX8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array3d_integer(invalue, outvalue, N, ngroup, dest)
    use mpi, only: MPI_INTEGER, MPI_SUM
    implicit none
    integer, intent(in)  :: invalue(:,:,:)
    integer, intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_INTEGER, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_INTEGER, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array3d_real8(invalue, outvalue, N, ngroup, dest)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array3d_complex8(invalue, outvalue, N, ngroup, dest)
    use mpi, only: MPI_COMPLEX8, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_COMPLEX8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_COMPLEX8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array4d_real8(invalue, outvalue, N, ngroup, dest)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_summation_array4d_complex8(invalue, outvalue, N, ngroup, dest)
    use mpi, only: MPI_COMPLEX8, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_COMPLEX8, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_COMPLEX8, MPI_SUM, ngroup, ierr))
    end if
  end subroutine

  subroutine comm_bcast_integer(val, ngroup, root)
    use mpi
    implicit none
    integer, intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, 1, MPI_INTEGER, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_real8(val, ngroup, root)
    use mpi
    implicit none
    real(8), intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, 1, MPI_REAL8, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_character(val, ngroup, root)
    use mpi
    implicit none
    character(*), intent(inout)        :: val
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, len(val), MPI_CHARACTER, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_logical(val, ngroup, root)
    use mpi
    implicit none
    logical, intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, 1, MPI_LOGICAL, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_array1d_integer(val, ngroup, root)
    use mpi
    implicit none
    integer, intent(inout)        :: val(:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_INTEGER, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_array1d_real8(val, ngroup, root)
    use mpi
    implicit none
    real(8), intent(inout)        :: val(:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_REAL8, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_array2d_integer(val, ngroup, root)
    use mpi
    implicit none
    integer, intent(inout)        :: val(:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_INTEGER, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_array2d_real8(val, ngroup, root)
    use mpi
    implicit none
    real(8), intent(inout)        :: val(:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_REAL8, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_array3d_real8(val, ngroup, root)
    use mpi
    implicit none
    real(8), intent(inout)        :: val(:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_REAL8, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_array3d_complex8(val, ngroup, root)
    use mpi
    implicit none
    complex(8), intent(inout)     :: val(:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_COMPLEX8, rank, ngroup, ierr))
  end subroutine

  subroutine comm_bcast_array1d_character(val, ngroup, root)
    use mpi
    implicit none
    character(*), intent(inout)        :: val(:)
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val)*len(val), MPI_CHARACTER, rank, ngroup, ierr))
  end subroutine

  subroutine comm_allgatherv_array1d_real8(invalue, outvalue, ncounts, displs, ngroup)
    use mpi, only: MPI_REAL8
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: ncounts(:)
    integer, intent(in)  :: displs(:)
    integer, intent(in)  :: ngroup
    integer :: ierr
    call MPI_Allgatherv(invalue,  size(invalue),   MPI_REAL8, &
                        outvalue, ncounts, displs, MPI_REAL8, &
                        ngroup, ierr)
    call error_check(ierr)
  end subroutine

  subroutine comm_get_min_array1d_real8(invalue, outvalue, N, ngroup)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_MIN, ngroup, ierr))
  end subroutine

  subroutine comm_get_max_array1d_real8(invalue, outvalue, N, ngroup)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_MAX, ngroup, ierr))
  end subroutine

  subroutine comm_get_maxloc(invalue, outvalue, ngroup)
    use mpi
    type(comm_maxloc_type), intent(in)  :: invalue
    type(comm_maxloc_type), intent(out) :: outvalue
    integer, intent(in)                 :: ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_DOUBLE_INT, MPI_MAXLOC, ngroup, ierr))
  end subroutine

  subroutine comm_logical_and_scalar(invalue, outvalue, ngroup)
    use mpi
    implicit none
    logical, intent(in)  :: invalue
    logical, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_LOGICAL, MPI_LAND, ngroup, ierr))
  end subroutine


  subroutine get_rank(comm, npid, nprocs)
    use mpi
    implicit none
    integer, intent(in)  :: comm
    integer, intent(out) :: npid, nprocs
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Comm_rank(comm,   npid, ierr))
    MPI_ERROR_CHECK(call MPI_Comm_size(comm, nprocs, ierr))
  end subroutine

  subroutine error_check(errcode)
    use mpi
    implicit none
    integer, intent(in) :: errcode
    character(MPI_MAX_ERROR_STRING) :: errstr
    integer                         :: retlen, ierr
    if (errcode /= MPI_SUCCESS) then
      call MPI_Error_string(errcode, errstr, retlen, ierr)
      print *, 'MPI Error:', errstr
      call comm_finalize
    end if
  end subroutine
end module
