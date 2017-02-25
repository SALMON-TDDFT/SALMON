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
module communication
  implicit none

  integer, public, parameter :: ROOT_PROCID = 0

  integer, public :: procid(2)
  integer, public :: nprocs(2)
  integer, public :: proc_group(2)

  public :: comm_init
  public :: comm_finalize
  public :: comm_is_root
  public :: comm_set_level2_group
  public :: comm_sync_all
  public :: comm_summation
  public :: comm_bcast
  public :: comm_get_min
  public :: comm_get_max

  type, public :: comm_maxloc_type
    real(8) :: val
    integer :: rank
  end type

  interface comm_summation
    ! scalar
    module procedure comm_summation_real8

    ! 1-D array
    module procedure comm_summation_array1d_real8

    ! 2-D array
    module procedure comm_summation_array2d_real8

    ! 3-D array
    module procedure comm_summation_array3d_real8
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

    ! 2-D array
    module procedure comm_bcast_array2d_integer
    module procedure comm_bcast_array2d_real8

    ! 3-D array
    module procedure comm_bcast_array3d_real8
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
    proc_group(:) = MPI_COMM_WORLD
    call get_rank(1)
    call get_rank(2)
  end subroutine

  subroutine comm_finalize
    use mpi
    implicit none
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Finalize(ierr))
  end subroutine

  function comm_is_root(level)
    implicit none
    integer, intent(in), optional :: level
    logical :: comm_is_root
    if (present(level)) then
      comm_is_root = procid(level) == ROOT_PROCID
    else
      comm_is_root = procid(1)     == ROOT_PROCID
    end if
  end function

  subroutine comm_set_level2_group(nprocs, key)
    use mpi
    implicit none
    integer, intent(in) :: nprocs, key
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Comm_split(proc_group(1), nprocs, key, proc_group(2), ierr))
    call get_rank(2)
  end subroutine

  subroutine comm_sync_all(level)
    use mpi
    implicit none
    integer, intent(in), optional :: level
    integer :: ierr
    if (present(level)) then
      MPI_ERROR_CHECK(call MPI_Barrier(level, ierr))
    else
      MPI_ERROR_CHECK(call MPI_Barrier(proc_group(1), ierr))
    end if
  end subroutine


  subroutine comm_summation_real8(invalue, outvalue, level)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue
    real(8), intent(out) :: outvalue
    integer, intent(in)  :: level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_REAL8, MPI_SUM, level, ierr))
  end subroutine

  subroutine comm_summation_array1d_real8(invalue, outvalue, N, level)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, level, ierr))
  end subroutine

  subroutine comm_summation_array2d_real8(invalue, outvalue, N, level)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:,:)
    real(8), intent(out) :: outvalue(:,:)
    integer, intent(in)  :: N, level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, level, ierr))
  end subroutine

  subroutine comm_summation_array3d_real8(invalue, outvalue, N, level)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:,:,:)
    real(8), intent(out) :: outvalue(:,:,:)
    integer, intent(in)  :: N, level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_SUM, level, ierr))
  end subroutine

  subroutine comm_bcast_integer(val, level, root)
    use mpi
    implicit none
    integer, intent(inout)        :: val
    integer, intent(in)           :: level
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, 1, MPI_INTEGER, rank, level, ierr))
  end subroutine

  subroutine comm_bcast_real8(val, level, root)
    use mpi
    implicit none
    real(8), intent(inout)        :: val
    integer, intent(in)           :: level
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, 1, MPI_REAL8, rank, level, ierr))
  end subroutine

  subroutine comm_bcast_character(val, level, root)
    use mpi
    implicit none
    character(*), intent(inout)        :: val
    integer,      intent(in)           :: level
    integer,      intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, len(val), MPI_CHARACTER, rank, level, ierr))
  end subroutine

  subroutine comm_bcast_logical(val, level, root)
    use mpi
    implicit none
    logical, intent(inout)        :: val
    integer, intent(in)           :: level
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, 1, MPI_LOGICAL, rank, level, ierr))
  end subroutine

  subroutine comm_bcast_array1d_integer(val, level, root)
    use mpi
    implicit none
    integer, intent(inout)        :: val(:)
    integer, intent(in)           :: level
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_INTEGER, rank, level, ierr))
  end subroutine

  subroutine comm_bcast_array1d_real8(val, level, root)
    use mpi
    implicit none
    real(8), intent(inout)        :: val(:)
    integer, intent(in)           :: level
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_REAL8, rank, level, ierr))
  end subroutine

  subroutine comm_bcast_array2d_integer(val, level, root)
    use mpi
    implicit none
    integer, intent(inout)        :: val(:,:)
    integer, intent(in)           :: level
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_INTEGER, rank, level, ierr))
  end subroutine

  subroutine comm_bcast_array2d_real8(val, level, root)
    use mpi
    implicit none
    real(8), intent(inout)        :: val(:,:)
    integer, intent(in)           :: level
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_REAL8, rank, level, ierr))
  end subroutine

  subroutine comm_bcast_array3d_real8(val, level, root)
    use mpi
    implicit none
    real(8), intent(inout)        :: val(:,:,:)
    integer, intent(in)           :: level
    integer, intent(in), optional :: root
    integer :: rank, ierr
    if (present(root)) then
      rank = root
    else
      rank = 0
    end if
    MPI_ERROR_CHECK(call MPI_Bcast(val, size(val), MPI_REAL8, rank, level, ierr))
  end subroutine

  subroutine comm_get_min_array1d_real8(invalue, outvalue, N, level)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_MIN, level, ierr))
  end subroutine

  subroutine comm_get_max_array1d_real8(invalue, outvalue, N, level)
    use mpi
    implicit none
    real(8), intent(in)  :: invalue(:)
    real(8), intent(out) :: outvalue(:)
    integer, intent(in)  :: N, level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, N, MPI_REAL8, MPI_MAX, level, ierr))
  end subroutine

  subroutine comm_get_maxloc(invalue, outvalue, level)
    use mpi
    type(comm_maxloc_type), intent(in)  :: invalue
    type(comm_maxloc_type), intent(out) :: outvalue
    integer, intent(in)                 :: level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_DOUBLE_INT, MPI_MAXLOC, level, ierr))
  end subroutine

  subroutine comm_logical_and_scalar(invalue, outvalue, level)
    use mpi
    implicit none
    logical, intent(in)  :: invalue
    logical, intent(out) :: outvalue
    integer, intent(in)  :: level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Allreduce(invalue, outvalue, 1, MPI_LOGICAL, MPI_LAND, level, ierr))
  end subroutine


  subroutine get_rank(level)
    use mpi
    implicit none
    integer, intent(in) :: level
    integer :: ierr
    MPI_ERROR_CHECK(call MPI_Comm_rank(proc_group(level), procid(level), ierr))
    MPI_ERROR_CHECK(call MPI_Comm_size(proc_group(level), nprocs(level), ierr))
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
