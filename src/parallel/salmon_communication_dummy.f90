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
!!! This is a dummy module of salmon_communication
module salmon_communication
  implicit none

  integer, public, parameter  :: ROOT_PROCID    = 0
  integer, public, parameter  :: COMM_PROC_NULL = int(z'0000BEEF')
  integer, private, parameter :: COMM_WORLD_ID  = int(z'000ABEEF')

  ! call once
  public :: comm_init
  public :: comm_finalize

  ! application stops when a following routines is called
  public :: comm_send
  public :: comm_recv
  public :: comm_exchange
  public :: comm_isend
  public :: comm_irecv
  public :: comm_wait
  public :: comm_wait_all
  public :: comm_send_init
  public :: comm_recv_init
  public :: comm_start_all

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

  private :: get_rank, abort_show_message

#define ABORT_MESSAGE(target,msg) if(target/=COMM_PROC_NULL) call abort_show_message(msg)
#define UNUSED_VARIABLE(VAR)      if(.false.) call salmon_unusedvar(VAR)

contains
  subroutine comm_init
    implicit none
    ! do nothing
  end subroutine

  subroutine comm_finalize
    implicit none
    ! do nothing
  end subroutine

  subroutine comm_get_globalinfo(ngid, npid, nprocs)
    implicit none
    integer, intent(out) :: ngid, npid, nprocs
    ngid = COMM_WORLD_ID
    call get_rank(ngid, npid, nprocs)
  end subroutine

  subroutine comm_get_groupinfo(ngid, npid, nprocs)
    implicit none
    integer, intent(in)  :: ngid
    integer, intent(out) :: npid, nprocs
    call get_rank(ngid, npid, nprocs)
  end subroutine

  function comm_create_group(ngid, nprocs, key) result(ngid_dst)
    implicit none
    integer, intent(in) :: ngid, nprocs, key
    integer :: ngid_dst
    ngid_dst = ngid + key * nprocs
  end function

  function comm_is_root(npid)
    implicit none
    integer, intent(in) :: npid
    logical :: comm_is_root
    comm_is_root = npid == ROOT_PROCID
  end function

  subroutine comm_sync_all(ngid)
    implicit none
    integer, intent(in), optional :: ngid
    UNUSED_VARIABLE(ngid)
    ! do nothing
  end subroutine


  subroutine comm_send_array5d_double(invalue, ndest, ntag, ngroup)
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_array5d_double")
  end subroutine

  subroutine comm_send_array5d_dcomplex(invalue, ndest, ntag, ngroup)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(ndest,"comm_send_array5d_dcomplex")
  end subroutine

  subroutine comm_recv_array5d_double(outvalue, nsrc, ntag, ngroup)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_array5d_double")
  end subroutine

  subroutine comm_recv_array5d_dcomplex(outvalue, nsrc, ntag, ngroup)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    ABORT_MESSAGE(nsrc,"comm_recv_array5d_dcomplex")
  end subroutine


  subroutine comm_exchange_array3d_double(invalue, ndest, outvalue, nsrc, ntag, ngroup)
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
  end subroutine

  subroutine comm_exchange_array3d_dcomplex(invalue, ndest, outvalue, nsrc, ntag, ngroup)
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
  end subroutine

  subroutine comm_exchange_array5d_double(invalue, ndest, outvalue, nsrc, ntag, ngroup)
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
  end subroutine

  subroutine comm_exchange_array5d_dcomplex(invalue, ndest, outvalue, nsrc, ntag, ngroup)
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
  end subroutine


  function comm_isend_array3d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
    ABORT_MESSAGE(ndest,"comm_isend_array3d_double")
  end function

  function comm_isend_array3d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
    ABORT_MESSAGE(ndest,"comm_isend_array3d_dcomplex")
  end function

  function comm_isend_array5d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
    ABORT_MESSAGE(ndest,"comm_isend_array5d_double")
  end function

  function comm_isend_array5d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
    ABORT_MESSAGE(ndest,"comm_isend_array5d_dcomplex")
  end function

  function comm_irecv_array3d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
    ABORT_MESSAGE(nsrc,"comm_irecv_array3d_double")
  end function

  function comm_irecv_array3d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
    ABORT_MESSAGE(nsrc,"comm_irecv_array3d_dcomplex")
  end function

  function comm_irecv_array5d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
    ABORT_MESSAGE(nsrc,"comm_irecv_array5d_double")
  end function

  function comm_irecv_array5d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
    ABORT_MESSAGE(nsrc,"comm_irecv_array5d_dcomplex")
  end function


  subroutine comm_wait(req)
    implicit none
    integer, intent(in) :: req
    UNUSED_VARIABLE(req)
    ! do nothing
  end subroutine

  subroutine comm_wait_all(reqs)
    implicit none
    integer, intent(in) :: reqs(:)
    UNUSED_VARIABLE(reqs)
    ! do nothing
  end subroutine


  function comm_send_init_array3d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ndest)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
  end function

  function comm_send_init_array3d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ndest)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
  end function

  function comm_send_init_array5d_double(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    real(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in) :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ndest)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
  end function

  function comm_send_init_array5d_dcomplex(invalue, ndest, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(in) :: invalue(:,:,:,:,:)
    integer, intent(in)    :: ndest, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(invalue)
    UNUSED_VARIABLE(ndest)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
  end function

  function comm_recv_init_array3d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
  end function

  function comm_recv_init_array3d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
  end function

  function comm_recv_init_array5d_double(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
  end function

  function comm_recv_init_array5d_dcomplex(outvalue, nsrc, ntag, ngroup) result(req)
    implicit none
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: nsrc, ntag, ngroup
    integer :: req
    UNUSED_VARIABLE(outvalue)
    UNUSED_VARIABLE(nsrc)
    UNUSED_VARIABLE(ntag)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(req)
  end function

  subroutine comm_start_all(reqs)
    implicit none
    integer, intent(in) :: reqs(:)
    UNUSED_VARIABLE(reqs)
  end subroutine


  subroutine comm_summation_integer(invalue, outvalue, ngroup, dest)
    implicit none
    integer, intent(in)  :: invalue
    integer, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_double(invalue, outvalue, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue
    real(8), intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_dcomplex(invalue, outvalue, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue
    complex(8), intent(out) :: outvalue
    integer, intent(in)     :: ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array1d_integer(invalue, outvalue, N, ngroup, dest)
    implicit none
    integer, intent(in)  :: invalue(:)
    integer, intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array1d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array1d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:)
    complex(8), intent(out) :: outvalue(:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array2d_integer(invalue, outvalue, N, ngroup, dest)
    implicit none
    integer, intent(in)  :: invalue(:,:)
    integer, intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array2d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:)
    real(8), intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array2d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:)
    complex(8), intent(out) :: outvalue(:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array3d_integer(invalue, outvalue, N, ngroup, dest)
    implicit none
    integer, intent(in)  :: invalue(:,:,:)
    integer, intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array3d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array3d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array4d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array4d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array5d_double(invalue, outvalue, N, ngroup, dest)
    implicit none
    real(8), intent(in)  :: invalue(:,:,:,:,:)
    real(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)  :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_summation_array5d_dcomplex(invalue, outvalue, N, ngroup, dest)
    implicit none
    complex(8), intent(in)  :: invalue(:,:,:,:,:)
    complex(8), intent(out) :: outvalue(:,:,:,:,:)
    integer, intent(in)     :: N, ngroup
    integer, optional, intent(in) :: dest
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(dest)
    outvalue = invalue
  end subroutine

  subroutine comm_bcast_integer(val, ngroup, root)
    implicit none
    integer, intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_character(val, ngroup, root)
    implicit none
    character(*), intent(inout)        :: val
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_logical(val, ngroup, root)
    implicit none
    logical, intent(inout)        :: val
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_array1d_integer(val, ngroup, root)
    implicit none
    integer, intent(inout)        :: val(:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_array1d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_array2d_integer(val, ngroup, root)
    implicit none
    integer, intent(inout)        :: val(:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_array2d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_array3d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine

  subroutine comm_bcast_array3d_dcomplex(val, ngroup, root)
    implicit none
    complex(8), intent(inout)     :: val(:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine
  
  subroutine comm_bcast_array4d_double(val, ngroup, root)
    implicit none
    real(8), intent(inout)        :: val(:,:,:,:)
    integer, intent(in)           :: ngroup
    integer, intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine
  
  subroutine comm_bcast_array1d_character(val, ngroup, root)
    implicit none
    character(*), intent(inout)        :: val(:)
    integer,      intent(in)           :: ngroup
    integer,      intent(in), optional :: root
    UNUSED_VARIABLE(val)
    UNUSED_VARIABLE(ngroup)
    UNUSED_VARIABLE(root)
    ! do nothing
  end subroutine


  subroutine comm_allgatherv_array1d_double(invalue, outvalue, ncounts, displs, ngroup)
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: ncounts(:)
    integer, intent(in)  :: displs(:)
    integer, intent(in)  :: ngroup
    UNUSED_VARIABLE(ngroup)
    outvalue(displs(1)+1:displs(1)+ncounts(1)) = invalue(1:ncounts(1)-1)
  end subroutine


  subroutine comm_get_min_array1d_double(invalue, outvalue, N, ngroup)
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    outvalue = invalue
  end subroutine

  subroutine comm_get_max_array1d_double(invalue, outvalue, N, ngroup)
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, ngroup
    UNUSED_VARIABLE(N)
    UNUSED_VARIABLE(ngroup)
    outvalue = invalue
  end subroutine

  subroutine comm_get_maxloc(invalue, outvalue, ngroup)
    implicit none
    type(comm_maxloc_type), intent(in)  :: invalue
    type(comm_maxloc_type), intent(out) :: outvalue
    integer, intent(in)                 :: ngroup
    UNUSED_VARIABLE(ngroup)
    outvalue = invalue
  end subroutine

  subroutine comm_logical_and_scalar(invalue, outvalue, ngroup)
    implicit none
    logical, intent(in)  :: invalue
    logical, intent(out) :: outvalue
    integer, intent(in)  :: ngroup
    UNUSED_VARIABLE(ngroup)
    outvalue = invalue
  end subroutine


  subroutine get_rank(comm, npid, nprocs)
    implicit none
    integer, intent(in)  :: comm
    integer, intent(out) :: npid, nprocs
    UNUSED_VARIABLE(comm)
    npid   = 0
    nprocs = 1
  end subroutine

  subroutine abort_show_message(msg)
    implicit none
    character(*), intent(in) :: msg
    print '(A,A)', msg, ': this subroutine must not called (it takes MPI)'
    stop
  end subroutine
end module
