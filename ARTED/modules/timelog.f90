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
module timelog
  use misc_routines, only: get_wtime
  implicit none

  integer,public,parameter :: LOG_DT_EVOLVE    = 0
  integer,public,parameter :: LOG_HPSI         = 1
  integer,public,parameter :: LOG_PSI_RHO      = 2
  integer,public,parameter :: LOG_HARTREE      = 3
  integer,public,parameter :: LOG_EXC_COR      = 4
  integer,public,parameter :: LOG_CURRENT      = 5
  integer,public,parameter :: LOG_TOTAL_ENERGY = 6
  integer,public,parameter :: LOG_ION_FORCE    = 7
  integer,public,parameter :: LOG_DT_EVOLVE_AC = 8
  integer,public,parameter :: LOG_K_SHIFT_WF   = 9
  integer,public,parameter :: LOG_OTHER        = 10
  integer,public,parameter :: LOG_ALLREDUCE    = 11

  integer,public,parameter :: LOG_CG           = 12
  integer,public,parameter :: LOG_DIAG         = 13
  integer,public,parameter :: LOG_SP_ENERGY    = 14
  integer,public,parameter :: LOG_GRAM_SCHMIDT = 15

  integer,public,parameter :: LOG_HPSI_INIT    = 16
  integer,public,parameter :: LOG_HPSI_STENCIL = 17
  integer,public,parameter :: LOG_HPSI_PSEUDO  = 18
  integer,public,parameter :: LOG_HPSI_UPDATE  = 19

  integer,public,parameter :: LOG_DYNAMICS     = 20

  public :: timelog_initialize
  public :: timelog_set, timelog_reset
  public :: timelog_get, timelog_thread_get
  public :: timelog_begin, timelog_end
  public :: timelog_thread_begin, timelog_thread_end
  public :: timelog_show_hour, timelog_show_min

  public :: timelog_reentrance_read, timelog_reentrance_write
  public :: timelog_write, timelog_thread_write


  integer,private,parameter   :: LOG_SIZE = 30
  real(8),private,allocatable :: log_time(:)
  real(8),private,allocatable :: log_temp(:)

  real(8),private,allocatable :: log_time_t(:,:)
  real(8),private,allocatable :: log_temp_t(:,:)

private
contains
  subroutine timelog_initialize
    use omp_lib, only: omp_get_max_threads
    implicit none
    allocate(log_time(0:LOG_SIZE - 1))
    allocate(log_temp(0:LOG_SIZE - 1))
    allocate(log_time_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    allocate(log_temp_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    call timelog_reset
  end subroutine

  subroutine timelog_set(e,t)
    implicit none
    integer,intent(in) :: e
    real(8),intent(in) :: t
    log_time(e) = t
    log_temp(e) = 0.d0
  end subroutine

  subroutine timelog_reset(e)
    implicit none
    integer,intent(in),optional :: e
    if(present(e)) then
      log_time  (e)   = 0.d0
      log_temp  (e)   = 0.d0
      log_time_t(e,:) = 0.d0
      log_temp_t(e,:) = 0.d0
    else
      log_time  (:)   = 0.d0
      log_temp  (:)   = 0.d0
      log_time_t(:,:) = 0.d0
      log_time_t(:,:) = 0.d0
    end if
  end subroutine

  subroutine timelog_reentrance_read(fd)
    implicit none
    integer,intent(in) :: fd
    read(fd) log_time(0:LOG_SIZE - 1)
    read(fd) log_temp(0:LOG_SIZE - 1)
  end subroutine

  subroutine timelog_reentrance_write(fd)
    implicit none
    integer,intent(in) :: fd
    write(fd) log_time(0:LOG_SIZE - 1)
    write(fd) log_temp(0:LOG_SIZE - 1)
  end subroutine

  subroutine timelog_begin(id)
    implicit none
    integer,intent(in) :: id
    log_temp(id) = get_wtime()
  end subroutine

  subroutine timelog_end(id)
    implicit none
    integer,intent(in) :: id
    log_time(id) = log_time(id) + get_wtime() - log_temp(id)
  end subroutine

  subroutine timelog_thread_begin(id)
    use omp_lib
    implicit none
    integer,intent(in) :: id
    integer :: tid
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timelog_begin(id)
    end if
    log_temp_t(id,tid) = get_wtime()
  end subroutine

  subroutine timelog_thread_end(id)
    use omp_lib
    implicit none
    integer,intent(in) :: id
    integer :: tid
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timelog_end(id)
    end if
    log_time_t(id,tid) = log_time_t(id,tid) + get_wtime() - log_temp_t(id,tid)
  end subroutine

  subroutine timelog_show_hour(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,hour
    time = log_time(id)
    hour = time / 3600
    write(*,*) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timelog_show_min(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,mini
    time = log_time(id)
    mini = time / 60
    write(*,*) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timelog_write(fd,str,id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,id
    real(8) :: time
    time = log_time(id)
    write(fd,*) str,time,'sec'
  end subroutine

  subroutine timelog_thread_write(fd,str,id)
    use omp_lib
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,id
    real(8) :: time
    integer :: i
    write(fd,*) str
    do i=0,omp_get_max_threads()-1
      time = log_time_t(id,i)
      write(fd,*) 'tid =',i,': ',time,'sec'
    end do
  end subroutine

  function timelog_get(id)
    implicit none
    integer,intent(in) :: id
    real(8)            :: timelog_get
    timelog_get = log_time(id)
  end function

  function timelog_thread_get(id,tid)
    implicit none
    integer,intent(in) :: id,tid
    real(8)            :: timelog_thread_get
    timelog_thread_get = log_time_t(id,tid)
  end function
end module timelog
