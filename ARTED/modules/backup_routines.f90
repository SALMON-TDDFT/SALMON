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
module backup_routines
  ! NOTE: `allocatable dummy array` feature requires Fortran2003 compiler.
  !       A feature is defined at `ISO TR 15581 Allocatable Enhancements`.
  implicit none

  public :: save_value
  public :: load_value
  public :: backup_value

  interface save_value
    module procedure save_character
    module procedure save_logical

    ! scalar
    module procedure save_real8
    module procedure save_complex8
    module procedure save_integer

    ! 1-D array
    module procedure save_array1d_real8
    module procedure save_array1d_complex8
    module procedure save_array1d_integer

    ! 2-D array
    module procedure save_array2d_real8
    module procedure save_array2d_complex8
    module procedure save_array2d_integer

    ! 3-D array
    module procedure save_array3d_real8
    module procedure save_array3d_complex8
    module procedure save_array3d_integer

    ! 4-D array
    module procedure save_array4d_real8
    module procedure save_array4d_complex8
    module procedure save_array4d_integer
  end interface

  interface load_value
    module procedure load_character
    module procedure load_logical

    ! scalar
    module procedure load_real8
    module procedure load_complex8
    module procedure load_integer

    ! 1-D array
    module procedure load_array1d_real8
    module procedure load_array1d_complex8
    module procedure load_array1d_integer

    ! 2-D array
    module procedure load_array2d_real8
    module procedure load_array2d_complex8
    module procedure load_array2d_integer

    ! 3-D array
    module procedure load_array3d_real8
    module procedure load_array3d_complex8
    module procedure load_array3d_integer

    ! 4-D array
    module procedure load_array4d_real8
    module procedure load_array4d_complex8
    module procedure load_array4d_integer
  end interface

  ! backup or restore
  interface backup_value
    module procedure backup_character
    module procedure backup_logical

    ! scalar
    module procedure backup_real8
    module procedure backup_complex8
    module procedure backup_integer

    ! 1-D array
    module procedure backup_array1d_real8
    module procedure backup_array1d_complex8
    module procedure backup_array1d_integer

    ! 2-D array
    module procedure backup_array2d_real8
    module procedure backup_array2d_complex8
    module procedure backup_array2d_integer

    ! 3-D array
    module procedure backup_array3d_real8
    module procedure backup_array3d_complex8
    module procedure backup_array3d_integer

    ! 4-D array
    module procedure backup_array4d_real8
    module procedure backup_array4d_complex8
    module procedure backup_array4d_integer
  end interface

contains
  subroutine save_character(iounit, val)
    implicit none
    integer, intent(in)      :: iounit
    character(*), intent(in) :: val
    write(iounit) val
  end subroutine

  subroutine load_character(iounit, val)
    implicit none
    integer, intent(in)       :: iounit
    character(*), intent(out) :: val
    read(iounit) val
  end subroutine

  subroutine backup_character(is_backup, iounit, val)
    implicit none
    logical, intent(in)         :: is_backup
    integer, intent(in)         :: iounit
    character(*), intent(inout) :: val
    if (is_backup) then
      call save_character(iounit, val)
    else
      call load_character(iounit, val)
    end if
  end subroutine

  subroutine save_logical(iounit, val)
    implicit none
    integer, intent(in) :: iounit
    logical, intent(in) :: val
    write(iounit) val
  end subroutine

  subroutine load_logical(iounit, val)
    implicit none
    integer, intent(in)  :: iounit
    logical, intent(out) :: val
    read(iounit) val
  end subroutine

  subroutine backup_logical(is_backup, iounit, val)
    implicit none
    logical, intent(in)    :: is_backup
    integer, intent(in)    :: iounit
    logical, intent(inout) :: val
    if (is_backup) then
      call save_logical(iounit, val)
    else
      call load_logical(iounit, val)
    end if
  end subroutine

  subroutine save_real8(iounit, val)
    implicit none
    integer, intent(in) :: iounit
    real(8), intent(in) :: val
    write(iounit) val
  end subroutine

  subroutine load_real8(iounit, val)
    implicit none
    integer, intent(in)  :: iounit
    real(8), intent(out) :: val
    read(iounit) val
  end subroutine

  subroutine backup_real8(is_backup, iounit, val)
    implicit none
    logical, intent(in)    :: is_backup
    integer, intent(in)    :: iounit
    real(8), intent(inout) :: val
    if (is_backup) then
      call save_real8(iounit, val)
    else
      call load_real8(iounit, val)
    end if
  end subroutine

  subroutine save_array1d_real8(iounit, val)
    implicit none
    integer, intent(in)              :: iounit
    real(8), allocatable, intent(in) :: val(:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val), ubound(val)
      write(iounit) val(:)
    end if
  end subroutine

  subroutine load_array1d_real8(iounit, val)
    implicit none
    integer, intent(in)               :: iounit
    real(8), allocatable, intent(out) :: val(:)
    logical :: f
    integer :: lb1, ub1
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      allocate(val(lb1:ub1))
      read(iounit) val(:)
    end if
  end subroutine

  subroutine backup_array1d_real8(is_backup, iounit, val)
    implicit none
    logical, intent(in)                 :: is_backup
    integer, intent(in)                 :: iounit
    real(8), allocatable, intent(inout) :: val(:)
    if (is_backup) then
      call save_array1d_real8(iounit, val)
    else
      call load_array1d_real8(iounit, val)
    end if
  end subroutine

  subroutine save_array2d_real8(iounit, val)
    implicit none
    integer, intent(in)              :: iounit
    real(8), allocatable, intent(in) :: val(:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) val(:,:)
    end if
  end subroutine

  subroutine load_array2d_real8(iounit, val)
    implicit none
    integer, intent(in)               :: iounit
    real(8), allocatable, intent(out) :: val(:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      allocate(val(lb1:ub1,lb2:ub2))
      read(iounit) val(:,:)
    end if
  end subroutine

  subroutine backup_array2d_real8(is_backup, iounit, val)
    implicit none
    logical, intent(in)                 :: is_backup
    integer, intent(in)                 :: iounit
    real(8), allocatable, intent(inout) :: val(:,:)
    if (is_backup) then
      call save_array2d_real8(iounit, val)
    else
      call load_array2d_real8(iounit, val)
    end if
  end subroutine

  subroutine save_array3d_real8(iounit, val)
    implicit none
    integer, intent(in)              :: iounit
    real(8), allocatable, intent(in) :: val(:,:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) lbound(val,3), ubound(val,3)
      write(iounit) val(:,:,:)
    end if
  end subroutine

  subroutine load_array3d_real8(iounit, val)
    implicit none
    integer, intent(in)               :: iounit
    real(8), allocatable, intent(out) :: val(:,:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    integer :: lb3, ub3
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      read(iounit) lb3, ub3
      allocate(val(lb1:ub1,lb2:ub2,lb3:ub3))
      read(iounit) val(:,:,:)
    end if
  end subroutine

  subroutine backup_array3d_real8(is_backup, iounit, val)
    implicit none
    logical, intent(in)                 :: is_backup
    integer, intent(in)                 :: iounit
    real(8), allocatable, intent(inout) :: val(:,:,:)
    if (is_backup) then
      call save_array3d_real8(iounit, val)
    else
      call load_array3d_real8(iounit, val)
    end if
  end subroutine

  subroutine save_array4d_real8(iounit, val)
    implicit none
    integer, intent(in)              :: iounit
    real(8), allocatable, intent(in) :: val(:,:,:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) lbound(val,3), ubound(val,3)
      write(iounit) lbound(val,4), ubound(val,4)
      write(iounit) val(:,:,:,:)
    end if
  end subroutine

  subroutine load_array4d_real8(iounit, val)
    implicit none
    integer, intent(in)               :: iounit
    real(8), allocatable, intent(out) :: val(:,:,:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    integer :: lb3, ub3
    integer :: lb4, ub4
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      read(iounit) lb3, ub3
      read(iounit) lb4, ub4
      allocate(val(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4))
      read(iounit) val(:,:,:,:)
    end if
  end subroutine

  subroutine backup_array4d_real8(is_backup, iounit, val)
    implicit none
    logical, intent(in)                 :: is_backup
    integer, intent(in)                 :: iounit
    real(8), allocatable, intent(inout) :: val(:,:,:,:)
    if (is_backup) then
      call save_array4d_real8(iounit, val)
    else
      call load_array4d_real8(iounit, val)
    end if
  end subroutine

  subroutine save_complex8(iounit, val)
    implicit none
    integer, intent(in)    :: iounit
    complex(8), intent(in) :: val
    write(iounit) val
  end subroutine

  subroutine load_complex8(iounit, val)
    implicit none
    integer, intent(in)     :: iounit
    complex(8), intent(out) :: val
    read(iounit) val
  end subroutine

  subroutine backup_complex8(is_backup, iounit, val)
    implicit none
    logical, intent(in)       :: is_backup
    integer, intent(in)       :: iounit
    complex(8), intent(inout) :: val
    if (is_backup) then
      call save_complex8(iounit, val)
    else
      call load_complex8(iounit, val)
    end if
  end subroutine

  subroutine save_array1d_complex8(iounit, val)
    implicit none
    integer, intent(in)                 :: iounit
    complex(8), allocatable, intent(in) :: val(:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val), ubound(val)
      write(iounit) val(:)
    end if
  end subroutine

  subroutine load_array1d_complex8(iounit, val)
    implicit none
    integer, intent(in)                  :: iounit
    complex(8), allocatable, intent(out) :: val(:)
    logical :: f
    integer :: lb1, ub1
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      allocate(val(lb1:ub1))
      read(iounit) val(:)
    end if
  end subroutine

  subroutine backup_array1d_complex8(is_backup, iounit, val)
    implicit none
    logical, intent(in)                    :: is_backup
    integer, intent(in)                    :: iounit
    complex(8), allocatable, intent(inout) :: val(:)
    if (is_backup) then
      call save_array1d_complex8(iounit, val)
    else
      call load_array1d_complex8(iounit, val)
    end if
  end subroutine

  subroutine save_array2d_complex8(iounit, val)
    implicit none
    integer, intent(in)                 :: iounit
    complex(8), allocatable, intent(in) :: val(:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) val(:,:)
    end if
  end subroutine

  subroutine load_array2d_complex8(iounit, val)
    implicit none
    integer, intent(in)                  :: iounit
    complex(8), allocatable, intent(out) :: val(:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      allocate(val(lb1:ub1,lb2:ub2))
      read(iounit) val(:,:)
    end if
  end subroutine

  subroutine backup_array2d_complex8(is_backup, iounit, val)
    implicit none
    logical, intent(in)                    :: is_backup
    integer, intent(in)                    :: iounit
    complex(8), allocatable, intent(inout) :: val(:,:)
    if (is_backup) then
      call save_array2d_complex8(iounit, val)
    else
      call load_array2d_complex8(iounit, val)
    end if
  end subroutine

  subroutine save_array3d_complex8(iounit, val)
    implicit none
    integer, intent(in)                 :: iounit
    complex(8), allocatable, intent(in) :: val(:,:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) lbound(val,3), ubound(val,3)
      write(iounit) val(:,:,:)
    end if
  end subroutine

  subroutine load_array3d_complex8(iounit, val)
    implicit none
    integer, intent(in)                  :: iounit
    complex(8), allocatable, intent(out) :: val(:,:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    integer :: lb3, ub3
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      read(iounit) lb3, ub3
      allocate(val(lb1:ub1,lb2:ub2,lb3:ub3))
      read(iounit) val(:,:,:)
    end if
  end subroutine

  subroutine backup_array3d_complex8(is_backup, iounit, val)
    implicit none
    logical, intent(in)                    :: is_backup
    integer, intent(in)                    :: iounit
    complex(8), allocatable, intent(inout) :: val(:,:,:)
    if (is_backup) then
      call save_array3d_complex8(iounit, val)
    else
      call load_array3d_complex8(iounit, val)
    end if
  end subroutine

  subroutine save_array4d_complex8(iounit, val)
    implicit none
    integer, intent(in)                 :: iounit
    complex(8), allocatable, intent(in) :: val(:,:,:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) lbound(val,3), ubound(val,3)
      write(iounit) lbound(val,4), ubound(val,4)
      write(iounit) val(:,:,:,:)
    end if
  end subroutine

  subroutine load_array4d_complex8(iounit, val)
    implicit none
    integer, intent(in)                  :: iounit
    complex(8), allocatable, intent(out) :: val(:,:,:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    integer :: lb3, ub3
    integer :: lb4, ub4
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      read(iounit) lb3, ub3
      read(iounit) lb4, ub4
      allocate(val(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4))
      read(iounit) val(:,:,:,:)
    end if
  end subroutine

  subroutine backup_array4d_complex8(is_backup, iounit, val)
    implicit none
    logical, intent(in)                    :: is_backup
    integer, intent(in)                    :: iounit
    complex(8), allocatable, intent(inout) :: val(:,:,:,:)
    if (is_backup) then
      call save_array4d_complex8(iounit, val)
    else
      call load_array4d_complex8(iounit, val)
    end if
  end subroutine

  subroutine save_integer(iounit, val)
    implicit none
    integer, intent(in) :: iounit
    integer, intent(in) :: val
    write(iounit) val
  end subroutine

  subroutine load_integer(iounit, val)
    implicit none
    integer, intent(in)  :: iounit
    integer, intent(out) :: val
    read(iounit) val
  end subroutine

  subroutine backup_integer(is_backup, iounit, val)
    implicit none
    logical, intent(in)    :: is_backup
    integer, intent(in)    :: iounit
    integer, intent(inout) :: val
    if (is_backup) then
      call save_integer(iounit, val)
    else
      call load_integer(iounit, val)
    end if
  end subroutine

  subroutine save_array1d_integer(iounit, val)
    implicit none
    integer, intent(in)              :: iounit
    integer, allocatable, intent(in) :: val(:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val), ubound(val)
      write(iounit) val(:)
    end if
  end subroutine

  subroutine load_array1d_integer(iounit, val)
    implicit none
    integer, intent(in)               :: iounit
    integer, allocatable, intent(out) :: val(:)
    logical :: f
    integer :: lb1, ub1
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      allocate(val(lb1:ub1))
      read(iounit) val(:)
    end if
  end subroutine

  subroutine backup_array1d_integer(is_backup, iounit, val)
    implicit none
    logical, intent(in)                 :: is_backup
    integer, intent(in)                 :: iounit
    integer, allocatable, intent(inout) :: val(:)
    if (is_backup) then
      call save_array1d_integer(iounit, val)
    else
      call load_array1d_integer(iounit, val)
    end if
  end subroutine

  subroutine save_array2d_integer(iounit, val)
    implicit none
    integer, intent(in)              :: iounit
    integer, allocatable, intent(in) :: val(:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) val(:,:)
    end if
  end subroutine

  subroutine load_array2d_integer(iounit, val)
    implicit none
    integer, intent(in)               :: iounit
    integer, allocatable, intent(out) :: val(:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      allocate(val(lb1:ub1,lb2:ub2))
      read(iounit) val(:,:)
    end if
  end subroutine

  subroutine backup_array2d_integer(is_backup, iounit, val)
    implicit none
    logical, intent(in)                 :: is_backup
    integer, intent(in)                 :: iounit
    integer, allocatable, intent(inout) :: val(:,:)
    if (is_backup) then
      call save_array2d_integer(iounit, val)
    else
      call load_array2d_integer(iounit, val)
    end if
  end subroutine

  subroutine save_array3d_integer(iounit, val)
    implicit none
    integer, intent(in)              :: iounit
    integer, allocatable, intent(in) :: val(:,:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) lbound(val,3), ubound(val,3)
      write(iounit) val(:,:,:)
    end if
  end subroutine

  subroutine load_array3d_integer(iounit, val)
    implicit none
    integer, intent(in)               :: iounit
    integer, allocatable, intent(out) :: val(:,:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    integer :: lb3, ub3
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      read(iounit) lb3, ub3
      allocate(val(lb1:ub1,lb2:ub2,lb3:ub3))
      read(iounit) val(:,:,:)
    end if
  end subroutine

  subroutine backup_array3d_integer(is_backup, iounit, val)
    implicit none
    logical, intent(in)                 :: is_backup
    integer, intent(in)                 :: iounit
    integer, allocatable, intent(inout) :: val(:,:,:)
    if (is_backup) then
      call save_array3d_integer(iounit, val)
    else
      call load_array3d_integer(iounit, val)
    end if
  end subroutine

  subroutine save_array4d_integer(iounit, val)
    implicit none
    integer, intent(in)              :: iounit
    integer, allocatable, intent(in) :: val(:,:,:,:)
    logical :: f
    f = allocated(val)
    write(iounit) f
    if (f) then
      write(iounit) lbound(val,1), ubound(val,1)
      write(iounit) lbound(val,2), ubound(val,2)
      write(iounit) lbound(val,3), ubound(val,3)
      write(iounit) lbound(val,4), ubound(val,4)
      write(iounit) val(:,:,:,:)
    end if
  end subroutine

  subroutine load_array4d_integer(iounit, val)
    implicit none
    integer, intent(in)               :: iounit
    integer, allocatable, intent(out) :: val(:,:,:,:)
    logical :: f
    integer :: lb1, ub1
    integer :: lb2, ub2
    integer :: lb3, ub3
    integer :: lb4, ub4
    read(iounit) f
    if (f) then
      read(iounit) lb1, ub1
      read(iounit) lb2, ub2
      read(iounit) lb3, ub3
      read(iounit) lb4, ub4
      allocate(val(lb1:ub1,lb2:ub2,lb3:ub3,lb4:ub4))
      read(iounit) val(:,:,:,:)
    end if
  end subroutine

  subroutine backup_array4d_integer(is_backup, iounit, val)
    implicit none
    logical, intent(in)                 :: is_backup
    integer, intent(in)                 :: iounit
    integer, allocatable, intent(inout) :: val(:,:,:,:)
    if (is_backup) then
      call save_array4d_integer(iounit, val)
    else
      call load_array4d_integer(iounit, val)
    end if
  end subroutine
end module
