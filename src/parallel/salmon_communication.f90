!
!  Copyright 2017 SALMON developers
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

#ifdef SALMON_USE_MPI
  use mpi, only: MPI_PROC_NULL
#endif

  implicit none

#ifdef SALMON_USE_MPI
  integer, public, parameter :: ROOT_PROCID = 0
  integer, public, parameter :: COMM_PROC_NULL = MPI_PROC_NULL
#else
  integer, public, parameter  :: ROOT_PROCID    = 0
  integer, public, parameter  :: COMM_PROC_NULL = int(z'0000BEEF')
  integer, private, parameter :: COMM_WORLD_ID  = COMM_PROC_NULL
  integer, private, parameter :: DEAD_BEAF      = int(z'0000DEAD')
#endif

  ! call once
  public :: comm_init
  public :: comm_finalize

  ! p2p communication
  ! application stops when a following routines is called in no-mpi environment
  public :: comm_send
  public :: comm_recv
  public :: comm_exchange

  ! p2p immediate communication
  ! application stops when a following routines is called in no-mpi environment
  public :: comm_isend
  public :: comm_irecv
  public :: comm_wait
  public :: comm_wait_all

  ! p2p persistent communication
  ! application stops when a following routines is called in no-mpi environment
  public :: comm_send_init
  public :: comm_recv_init
  public :: comm_start_all

  ! collective communication
  public :: comm_sync_all
  public :: comm_summation
  public :: comm_bcast
  public :: comm_allgatherv ! not implemented in no-mpi environment
  public :: comm_alltoall
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
    module procedure comm_send_array5d_double
    module procedure comm_send_array5d_dcomplex
  end interface

  interface comm_recv
    ! 5-D array
    module procedure comm_recv_array5d_double
    module procedure comm_recv_array5d_dcomplex
  end interface

  interface comm_exchange
    ! 3-D array
    module procedure comm_exchange_array3d_double
    module procedure comm_exchange_array3d_dcomplex

    ! 5-D array
    module procedure comm_exchange_array5d_double
    module procedure comm_exchange_array5d_dcomplex
  end interface

  interface comm_isend
    ! 3-D array
    module procedure comm_isend_array3d_double
    module procedure comm_isend_array3d_dcomplex

    ! 5-D array
    module procedure comm_isend_array5d_double
    module procedure comm_isend_array5d_dcomplex
  end interface

  interface comm_irecv
    ! 3-D array
    module procedure comm_irecv_array3d_double
    module procedure comm_irecv_array3d_dcomplex

    ! 5-D array
    module procedure comm_irecv_array5d_double
    module procedure comm_irecv_array5d_dcomplex
  end interface

  interface comm_send_init
    ! 3-D array
    module procedure comm_send_init_array3d_double
    module procedure comm_send_init_array3d_dcomplex

    ! 5-D array
    module procedure comm_send_init_array5d_double
    module procedure comm_send_init_array5d_dcomplex
  end interface

  interface comm_recv_init
    ! 3-D array
    module procedure comm_recv_init_array3d_double
    module procedure comm_recv_init_array3d_dcomplex

    ! 5-D array
    module procedure comm_recv_init_array5d_double
    module procedure comm_recv_init_array5d_dcomplex
  end interface

  interface comm_summation
    ! scalar
    module procedure comm_summation_integer
    module procedure comm_summation_double
    module procedure comm_summation_dcomplex

    ! 1-D array
    module procedure comm_summation_array1d_integer
    module procedure comm_summation_array1d_double
    module procedure comm_summation_array1d_dcomplex

    ! 2-D array
    module procedure comm_summation_array2d_integer
    module procedure comm_summation_array2d_double
    module procedure comm_summation_array2d_dcomplex

    ! 3-D array
    module procedure comm_summation_array3d_integer
    module procedure comm_summation_array3d_double
    module procedure comm_summation_array3d_dcomplex

    ! 4-D array
    module procedure comm_summation_array4d_double
    module procedure comm_summation_array4d_dcomplex

    ! 5-D array
    module procedure comm_summation_array5d_double
    module procedure comm_summation_array5d_dcomplex
  end interface

  interface comm_bcast
    ! scalar
    module procedure comm_bcast_integer
    module procedure comm_bcast_double
    module procedure comm_bcast_character
    module procedure comm_bcast_logical

    ! 1-D array
    module procedure comm_bcast_array1d_integer
    module procedure comm_bcast_array1d_double
    module procedure comm_bcast_array1d_character

    ! 2-D array
    module procedure comm_bcast_array2d_integer
    module procedure comm_bcast_array2d_double

    ! 3-D array
    module procedure comm_bcast_array3d_double
    module procedure comm_bcast_array3d_dcomplex
  
    ! 4-D array
    module procedure comm_bcast_array4d_double
    ! module procedure comm_bcast_array3d_dcomplex
    !! TODO: create broadcast routine for rank-4 tensor later ...
  end interface

  interface comm_allgatherv
    ! 1-D array
    module procedure comm_allgatherv_array1d_double
  end interface

  interface comm_alltoall
    ! 1-D array
    module procedure comm_alltoall_array1d_complex
  end interface

  interface comm_get_min
    ! 1-D array
    module procedure comm_get_min_array1d_double
  end interface

  interface comm_get_max
    ! scalar
    module procedure comm_get_maxloc

    ! 1-D array
    module procedure comm_get_max_array1d_double
  end interface

  interface comm_logical_and
    ! scalar
    module procedure comm_logical_and_scalar
  end interface

  private :: get_rank, error_check, abort_show_message

#define MPI_ERROR_CHECK(x) x; call error_check(ierr)
#define ABORT_MESSAGE(target,msg) if(target/=COMM_PROC_NULL) call abort_show_message(msg)
#define UNUSED_VARIABLE(VAR)      if(.false.) call salmon_unusedvar(VAR)

contains
  subroutine comm_init
    implicit none
#ifdef SALMON_USE_MPI
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Init(ierr))
#endif
  end subroutine

  subroutine comm_finalize
    implicit none
#ifdef SALMON_USE_MPI
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Finalize(ierr))
#endif
  end subroutine

  subroutine comm_get_globalinfo(ngid, npid, nprocs)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_COMM_WORLD
    implicit none
    integer, intent(out) :: ngid, npid, nprocs
    ngid = MPI_COMM_WORLD
    call get_rank(ngid, npid, nprocs)
#else
    implicit none
    integer, intent(out) :: ngid, npid, nprocs
    ngid = COMM_WORLD_ID
    call get_rank(ngid, npid, nprocs)
#endif
  end subroutine

  subroutine comm_get_groupinfo(ngid, npid, nprocs)
#ifdef SALMON_USE_MPI
    implicit none
    integer, intent(in)  :: ngid
    integer, intent(out) :: npid, nprocs
    call get_rank(ngid, npid, nprocs)
#else
    implicit none
    integer, intent(in)  :: ngid
    integer, intent(out) :: npid, nprocs
    call get_rank(ngid, npid, nprocs)
#endif
  end subroutine

  function comm_create_group(ngid, nprocs, key) result(ngid_dst)
#ifdef SALMON_USE_MPI
    implicit none
    integer, intent(in) :: ngid, nprocs, key
    integer :: ngid_dst, ierr
    MPI_ERROR_CHECK(call MPI_Comm_split(ngid, nprocs, key, ngid_dst, ierr))
#else
    implicit none
    integer, intent(in) :: ngid, nprocs, key
    integer :: ngid_dst
    UNUSED_VARIABLE(key)
    UNUSED_VARIABLE(nprocs)
    ngid_dst = ngid
#endif
  end function

  function comm_is_root(npid)
#ifdef SALMON_USE_MPI
    implicit none
    integer, intent(in) :: npid
    logical :: comm_is_root
    comm_is_root = npid == ROOT_PROCID
#else
    implicit none
    integer, intent(in) :: npid
    logical :: comm_is_root
    comm_is_root = npid == ROOT_PROCID
#endif
  end function

  subroutine comm_sync_all(ngid)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_COMM_WORLD
    implicit none
    integer, intent(in), optional :: ngid
    integer :: ierr
    if (present(ngid)) then
      MPI_ERROR_CHECK(call MPI_Barrier(ngid, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Barrier(MPI_COMM_WORLD, ierr))
    end if
#else
    implicit none
    integer, intent(in), optional :: ngid
    UNUSED_VARIABLE(ngid)
    if (present(ngid)) then
      ABORT_MESSAGE(ngid,"comm_sync_all")
    else
      ! do nothing
    end if
#endif
  end subroutine


  subroutine comm_send_array5d_double(invalue, ndest, ntag, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Send(invalue, size(invalue), MPI_DOUBLE_PRECISION, ndest, ntag, ngroup, ierr))
#else
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_array5d_double")
#endif
  end subroutine

  subroutine comm_send_array5d_dcomplex(invalue, ndest, ntag, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_STATUS_SIZE
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Send(invalue, size(invalue), MPI_DOUBLE_COMPLEX, ndest, ntag, ngroup, ierr))
#else
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_array5d_dcomplex")
#endif
  end subroutine

  subroutine comm_recv_array5d_double(outvalue, nsrc, ntag, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    MPI_ERROR_CHECK(call MPI_Recv(outvalue, size(outvalue), MPI_DOUBLE_PRECISION, nsrc, ntag, ngroup, istatus, ierr))
#else
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_array5d_double")
#endif
  end subroutine

  subroutine comm_recv_array5d_dcomplex(outvalue, nsrc, ntag, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_STATUS_SIZE
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    MPI_ERROR_CHECK(call MPI_Recv(outvalue, size(outvalue), MPI_DOUBLE_COMPLEX, nsrc, ntag, ngroup, istatus, ierr))
#else
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_array5d_dcomplex")
#endif
  end subroutine


  subroutine comm_exchange_array3d_double(invalue, ndest, outvalue, nsrc, ntag, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ndest, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    call MPI_Sendrecv(invalue,  size(invalue),  MPI_DOUBLE_PRECISION, ndest, ntag, &
                      outvalue, size(outvalue), MPI_DOUBLE_PRECISION, nsrc,  ntag, &
                      ngroup, istatus, ierr)
    call error_check(ierr)
#else
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_exchange_array3d_double")
#endif
  end subroutine

  subroutine comm_exchange_array3d_dcomplex(invalue, ndest, outvalue, nsrc, ntag, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_STATUS_SIZE
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ndest, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    call MPI_Sendrecv(invalue,  size(invalue),  MPI_DOUBLE_COMPLEX, ndest, ntag, &
                      outvalue, size(outvalue), MPI_DOUBLE_COMPLEX, nsrc,  ntag, &
                      ngroup, istatus, ierr)
    call error_check(ierr)
#else
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_exchange_array3d_dcomplex")
#endif
  end subroutine

  subroutine comm_exchange_array5d_double(invalue, ndest, outvalue, nsrc, ntag, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_SIZE
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ndest, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    call MPI_Sendrecv(invalue,  size(invalue),  MPI_DOUBLE_PRECISION, ndest, ntag, &
                      outvalue, size(outvalue), MPI_DOUBLE_PRECISION, nsrc,  ntag, &
                      ngroup, istatus, ierr)
    call error_check(ierr)
#else
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_exchange_array5d_double")
#endif
  end subroutine

  subroutine comm_exchange_array5d_dcomplex(invalue, ndest, outvalue, nsrc, ntag, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_STATUS_SIZE
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ndest, ntag, ngroup
    integer :: ierr,istatus(MPI_STATUS_SIZE)
    call MPI_Sendrecv(invalue,  size(invalue),  MPI_DOUBLE_COMPLEX, ndest, ntag, &
                      outvalue, size(outvalue), MPI_DOUBLE_COMPLEX, nsrc,  ntag, &
                      ngroup, istatus, ierr)
    call error_check(ierr)
#else
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_exchange_array5d_dcomplex")
#endif
  end subroutine


  function comm_isend_array3d_double(invalue, ndest, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
    implicit none
    real(8), intent(in) :: invalue(:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: ierr, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Isend(invalue, size(invalue), MPI_DOUBLE_PRECISION, ndest, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    real(8), intent(in) :: invalue(:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_isend_array3d_double")
    req = 0
#endif
  end function

  function comm_isend_array3d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE
    implicit none
    complex(8), intent(in) :: invalue(:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: ierr, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Isend(invalue, size(invalue), MPI_DOUBLE_COMPLEX, ndest, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    complex(8), intent(in) :: invalue(:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_isend_array3d_dcomplex")
    req = 0
#endif
  end function

  function comm_isend_array5d_double(invalue, ndest, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: ierr, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Isend(invalue, size(invalue), MPI_DOUBLE_PRECISION, ndest, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_isend_array5d_double")
    req = 0
#endif
  end function

  function comm_isend_array5d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: ierr, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Isend(invalue, size(invalue), MPI_DOUBLE_COMPLEX, ndest, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_isend_array5d_dcomplex")
    req = 0
#endif
  end function

  function comm_irecv_array3d_double(outvalue, nsrc, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
    implicit none
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: ierr, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Irecv(outvalue, size(outvalue), MPI_DOUBLE_PRECISION, nsrc, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_irecv_array3d_double")
    req = 0
#endif
  end function

  function comm_irecv_array3d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: ierr, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Irecv(outvalue, size(outvalue), MPI_DOUBLE_COMPLEX, nsrc, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_irecv_array3d_dcomplex")
    req = 0
#endif
  end function

  function comm_irecv_array5d_double(outvalue, nsrc, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: ierr, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Irecv(outvalue, size(outvalue), MPI_DOUBLE_PRECISION, nsrc, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_irecv_array5d_double")
    req = 0
#endif
  end function

  function comm_irecv_array5d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: ierr, req
    logical :: flag
    MPI_ERROR_CHECK(call MPI_Irecv(outvalue, size(outvalue), MPI_DOUBLE_COMPLEX, nsrc, ntag, ngroup, req, ierr))
    MPI_ERROR_CHECK(call MPI_Test(req, flag, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_irecv_array5d_dcomplex")
    req = 0
#endif
  end function


  subroutine comm_wait(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_STATUS_IGNORE
    implicit none
    integer, intent(in) :: req
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Wait(req, MPI_STATUS_IGNORE, ierr))
#else
    implicit none
    integer, intent(in) :: req
    ABORT_MESSAGE(req,"comm_wait")
#endif
  end subroutine

  subroutine comm_wait_all(reqs)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_STATUSES_IGNORE
    implicit none
    integer, intent(in) :: reqs(:)
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Waitall(size(reqs), reqs, MPI_STATUSES_IGNORE, ierr))
#else
    implicit none
    integer, intent(in) :: reqs(:)
    UNUSED_VARIABLE(reqs)
    ! do nothing
#endif
  end subroutine

  function comm_send_init_array3d_double(invalue, ndest, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
    implicit none
    real(8), intent(in) :: invalue(:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: ierr, req
    MPI_ERROR_CHECK(call MPI_Send_init(invalue, size(invalue), MPI_DOUBLE_PRECISION, ndest, ntag, ngroup, req, ierr))
#else
    implicit none
    real(8), intent(in) :: invalue(:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array3d_double")
    req = DEAD_BEAF
#endif
  end function

  function comm_send_init_array3d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX
    implicit none
    complex(8), intent(in) :: invalue(:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: ierr, req
    MPI_ERROR_CHECK(call MPI_Send_init(invalue, size(invalue), MPI_DOUBLE_COMPLEX, ndest, ntag, ngroup, req, ierr))
#else
    implicit none
    complex(8), intent(in) :: invalue(:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array3d_dcomplex")
    req = DEAD_BEAF
#endif
  end function

  function comm_send_init_array5d_double(invalue, ndest, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: ierr, req
    MPI_ERROR_CHECK(call MPI_Send_init(invalue, size(invalue), MPI_DOUBLE_PRECISION, ndest, ntag, ngroup, req, ierr))
#else
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array5d_double")
    req = DEAD_BEAF
#endif
  end function

  function comm_send_init_array5d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: ierr, req
    MPI_ERROR_CHECK(call MPI_Send_init(invalue, size(invalue), MPI_DOUBLE_COMPLEX, ndest, ntag, ngroup, req, ierr))
#else
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_init_array5d_dcomplex")
    req = DEAD_BEAF
#endif
  end function

  function comm_recv_init_array3d_double(outvalue, nsrc, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
    implicit none
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: ierr, req
    MPI_ERROR_CHECK(call MPI_Recv_init(outvalue, size(outvalue), MPI_DOUBLE_PRECISION, nsrc, ntag, ngroup, req, ierr))
#else
    implicit none
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array3d_double")
    req = DEAD_BEAF
#endif
  end function

  function comm_recv_init_array3d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: ierr, req
    MPI_ERROR_CHECK(call MPI_Recv_init(outvalue, size(outvalue), MPI_DOUBLE_COMPLEX, nsrc, ntag, ngroup, req, ierr))
#else
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array3d_dcomplex")
    req = DEAD_BEAF
#endif
  end function

  function comm_recv_init_array5d_double(outvalue, nsrc, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: ierr, req
    MPI_ERROR_CHECK(call MPI_Recv_init(outvalue, size(outvalue), MPI_DOUBLE_PRECISION, nsrc, ntag, ngroup, req, ierr))
#else
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array5d_double")
    req = DEAD_BEAF
#endif
  end function

  function comm_recv_init_array5d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: ierr, req
    MPI_ERROR_CHECK(call MPI_Recv_init(outvalue, size(outvalue), MPI_DOUBLE_COMPLEX, nsrc, ntag, ngroup, req, ierr))
#else
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_init_array5d_dcomplex")
    req = DEAD_BEAF
#endif
  end function

  subroutine comm_start_all(reqs)
#ifdef SALMON_USE_MPI
    implicit none
    integer, intent(in) :: reqs(:)
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Startall(size(reqs), reqs, ierr))
#else
    implicit none
    integer, intent(in) :: reqs(:)
    UNUSED_VARIABLE(reqs)
    ! do nothing
#endif
  end subroutine


  subroutine comm_summation_integer(invalue, outvalue, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_INTEGER, MPI_SUM
    implicit none
    integer, intent(in)  :: invalue
    integer, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, 1, MPI_INTEGER, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_INTEGER, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    integer, intent(in)  :: invalue
    integer, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_integer")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_double(invalue, outvalue, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM
    implicit none
    real(8), intent(in)  :: invalue
    real(8), intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    real(8), intent(in)  :: invalue
    real(8), intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_double")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_dcomplex(invalue, outvalue, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue
    complex(8), intent(out) :: outvalue
    integer, intent(in)     :: ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    complex(8), intent(in)  :: invalue
    complex(8), intent(out) :: outvalue
    integer, intent(in)     :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_dcomplex")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array1d_integer(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
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
#else
    implicit none
    integer, intent(in)  :: invalue(:)
    integer, intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array1d_integer")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array1d_double(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array1d_double")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array1d_dcomplex(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:)
    complex(8), intent(out) :: outvalue(:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    complex(8), intent(in)  :: invalue(:)
    complex(8), intent(out) :: outvalue(:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array1d_dcomplex")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array2d_integer(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
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
#else
    implicit none
    integer, intent(in)  :: invalue(:,:)
    integer, intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array2d_integer")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array2d_double(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM
    implicit none
    real(8), intent(in)  :: invalue(:,:)
    real(8), intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    real(8), intent(in)  :: invalue(:,:)
    real(8), intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array2d_double")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array2d_dcomplex(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:,:)
    complex(8), intent(out) :: outvalue(:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    complex(8), intent(in)  :: invalue(:,:)
    complex(8), intent(out) :: outvalue(:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array2d_dcomplex")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array3d_integer(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
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
#else
    implicit none
    integer, intent(in)  :: invalue(:,:,:)
    integer, intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array3d_integer")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array3d_double(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array3d_double")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array3d_dcomplex(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array3d_dcomplex")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array4d_double(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array4d_double")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array4d_dcomplex(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array4d_dcomplex")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array5d_double(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array5d_double")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_summation_array5d_dcomplex(invalue, outvalue, N, ngroup, dest)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX, MPI_SUM
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    integer :: ierr
    if (present(dest)) then
      MPI_ERROR_CHECK(call MPI_Reduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, dest, ngroup, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_COMPLEX, MPI_SUM, ngroup, ierr))
    end if
#else
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(dest)
    ABORT_MESSAGE(ngroup,"comm_summation_array5d_dcomplex")
    outvalue = invalue
#endif
  end subroutine


  subroutine comm_bcast_integer(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_INTEGER
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
#else
    implicit none
    integer, intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_integer")
#endif
  end subroutine

  subroutine comm_bcast_double(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
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
    MPI_ERROR_CHECK(call MPI_Bcast(val, 1, MPI_DOUBLE_PRECISION, rank, ngroup, ierr))
#else
    implicit none
    real(8), intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_double")
#endif
  end subroutine

  subroutine comm_bcast_character(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_CHARACTER
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
#else
    implicit none
    character(*), intent(inout)        :: val
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_character")
#endif
  end subroutine

  subroutine comm_bcast_logical(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_LOGICAL
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
#else
    implicit none
    logical, intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_logical")
#endif
  end subroutine

  subroutine comm_bcast_array1d_integer(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_INTEGER
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
#else
    implicit none
    integer, intent(inout)        :: val(:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_array1d_integer")
#endif
  end subroutine

  subroutine comm_bcast_array1d_double(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
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
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_DOUBLE_PRECISION, rank, ngroup, ierr))
#else
    implicit none
    real(8), intent(inout)        :: val(:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_array1d_double")
#endif
  end subroutine

  subroutine comm_bcast_array2d_integer(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_INTEGER
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
#else
    implicit none
    integer, intent(inout)        :: val(:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_array2d_integer")
#endif
  end subroutine

  subroutine comm_bcast_array2d_double(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
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
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_DOUBLE_PRECISION, rank, ngroup, ierr))
#else
    implicit none
    real(8), intent(inout)        :: val(:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_array2d_double")
#endif
  end subroutine

  subroutine comm_bcast_array3d_double(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
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
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_DOUBLE_PRECISION, rank, ngroup, ierr))
#else
    implicit none
    real(8), intent(inout)        :: val(:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_array3d_double")
#endif
  end subroutine
  
  subroutine comm_bcast_array4d_double(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
    implicit none
    real(8), intent(inout)        :: val(:,:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_DOUBLE_PRECISION, rank, ngroup, ierr))
#else
    implicit none
    real(8), intent(inout)        :: val(:,:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_array4d_double")
#endif
  end subroutine

  subroutine comm_bcast_array3d_dcomplex(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX
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
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_DOUBLE_COMPLEX, rank, ngroup, ierr))
#else
    implicit none
    complex(8), intent(inout)     :: val(:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_array3d_dcomplex")
#endif
  end subroutine

  subroutine comm_bcast_array1d_character(val, ngroup, root)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_CHARACTER
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
#else
    implicit none
    character(*), intent(inout)        :: val(:)
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(root)
    ABORT_MESSAGE(ngroup,"comm_bcast_array1d_character")
#endif
  end subroutine


  subroutine comm_allgatherv_array1d_double(invalue, outvalue, ncounts, displs, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: ncounts(:)
    integer, intent(in)  :: displs(:)
    integer, intent(in)  :: ngroup
    integer :: ierr
    call MPI_Allgatherv(invalue,  size(invalue),   MPI_DOUBLE_PRECISION, &
                        outvalue, ncounts, displs, MPI_DOUBLE_PRECISION, &
                        ngroup, ierr)
    call error_check(ierr)
#else
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: ncounts(:)
    integer, intent(in)  :: displs(:)
    integer, intent(in)  :: ngroup
    ABORT_MESSAGE(ngroup,"comm_allgatherv_array1d_double")
    outvalue(displs(1)+1:displs(1)+ncounts(1)) = invalue(1:ncounts(1)-1)
#endif
  end subroutine


  subroutine comm_alltoall_array1d_complex(invalue, outvalue, ngroup, ncount)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_COMPLEX
    implicit none
    complex(8), intent(in)  :: invalue(:)
    complex(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: ngroup
    integer, intent(in)  :: ncount
    integer :: ierr
    call MPI_Alltoall(invalue,  ncount,          MPI_DOUBLE_COMPLEX, &
                      outvalue, ncount,          MPI_DOUBLE_COMPLEX, &
                      ngroup, ierr)
    call error_check(ierr)
#else
    implicit none
    complex(8), intent(in)  :: invalue(:)
    complex(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: ngroup
    integer, intent(in)  :: ncount
    UNUSED_VARIABLE(ncount)
    ABORT_MESSAGE(ngroup,"comm_alltoall_array1d_complex")
    outvalue = invalue
#endif
  end subroutine
 
 
  subroutine comm_get_min_array1d_double(invalue, outvalue, N, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_MIN
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_MIN, ngroup, ierr))
#else
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    UNUSED_VARIABLE(N)
    ABORT_MESSAGE(ngroup,"comm_get_min_array1d_double")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_get_max_array1d_double(invalue, outvalue, N, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_MAX
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_DOUBLE_PRECISION, MPI_MAX, ngroup, ierr))
#else
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    UNUSED_VARIABLE(N)
    ABORT_MESSAGE(ngroup,"comm_get_max_array1d_double")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_get_maxloc(invalue, outvalue, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_DOUBLE_INT, MPI_MAXLOC
    type(comm_maxloc_type), intent(in)  :: invalue
    type(comm_maxloc_type), intent(out) :: outvalue
    integer, intent(in)                 :: ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_DOUBLE_INT, MPI_MAXLOC, ngroup, ierr))
#else
    implicit none
    type(comm_maxloc_type), intent(in)  :: invalue
    type(comm_maxloc_type), intent(out) :: outvalue
    integer, intent(in)                 :: ngroup
    ABORT_MESSAGE(ngroup,"comm_get_maxloc")
    outvalue = invalue
#endif
  end subroutine

  subroutine comm_logical_and_scalar(invalue, outvalue, ngroup)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_LOGICAL, MPI_LAND
    implicit none
    logical, intent(in)  :: invalue
    logical, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_LOGICAL, MPI_LAND, ngroup, ierr))
#else
    implicit none
    logical, intent(in)  :: invalue
    logical, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    ABORT_MESSAGE(ngroup,"comm_logical_and_scalar")
    outvalue = invalue
#endif
  end subroutine


  subroutine get_rank(comm, npid, nprocs)
#ifdef SALMON_USE_MPI
    implicit none
    integer, intent(in)  :: comm
    integer, intent(out) :: npid, nprocs
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Comm_rank(comm,   npid, ierr))
    MPI_ERROR_CHECK(call MPI_Comm_size(comm, nprocs, ierr))
#else
    implicit none
    integer, intent(in)  :: comm
    integer, intent(out) :: npid, nprocs
    UNUSED_VARIABLE(comm)
    npid   = 0
    nprocs = 1
#endif
  end subroutine

  subroutine error_check(errcode)
#ifdef SALMON_USE_MPI
    use mpi, only: MPI_MAX_ERROR_STRING, MPI_SUCCESS
    implicit none
    integer, intent(in) :: errcode
    character(MPI_MAX_ERROR_STRING) :: errstr
    integer                         :: retlen, ierr
    if (errcode /= MPI_SUCCESS) then
      call MPI_Error_string(errcode, errstr, retlen, ierr)
      print *, 'MPI Error:', errstr
      call comm_finalize
    end if
#else
    implicit none
    integer, intent(in) :: errcode
    UNUSED_VARIABLE(errcode)
    call abort_show_message('error_check')
#endif
  end subroutine

  subroutine abort_show_message(msg)
    implicit none
    character(*), intent(in) :: msg
    print '(A,A)', msg, ': this subroutine must not called (it takes MPI)'
#ifdef __INTEL_COMPILER
    call tracebackqq
#endif
    stop
  end subroutine
  
end module
