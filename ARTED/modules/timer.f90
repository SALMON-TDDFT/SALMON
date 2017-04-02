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
module timer
  use misc_routines, only: get_wtime
  implicit none

  ! Calculation
  integer,public,parameter :: LOG_ALL          = 0
  integer,public,parameter :: LOG_STATIC       = 1
  integer,public,parameter :: LOG_GROUND_STATE = 2
  integer,public,parameter :: LOG_DYNAMICS     = 3
  integer,public,parameter :: LOG_IO           = 7

  ! GS routines
  integer,public,parameter :: LOG_CG           = 10
  integer,public,parameter :: LOG_DIAG         = 11
  integer,public,parameter :: LOG_SP_ENERGY    = 12
  integer,public,parameter :: LOG_GRAM_SCHMIDT = 13

  ! GS and RT routines
  integer,public,parameter :: LOG_DT_EVOLVE    = 20
  integer,public,parameter :: LOG_HPSI         = 21
  integer,public,parameter :: LOG_PSI_RHO      = 22
  integer,public,parameter :: LOG_HARTREE      = 23
  integer,public,parameter :: LOG_EXC_COR      = 24
  integer,public,parameter :: LOG_CURRENT      = 25
  integer,public,parameter :: LOG_TOTAL_ENERGY = 26
  integer,public,parameter :: LOG_ION_FORCE    = 27
  integer,public,parameter :: LOG_DT_EVOLVE_AC = 28
  integer,public,parameter :: LOG_K_SHIFT_WF   = 29
  integer,public,parameter :: LOG_OTHER        = 30

  ! Hamiltonian
  integer,public,parameter :: LOG_HPSI_INIT    = 35
  integer,public,parameter :: LOG_HPSI_STENCIL = 36
  integer,public,parameter :: LOG_HPSI_PSEUDO  = 37
  integer,public,parameter :: LOG_HPSI_UPDATE  = 38

  ! Communication
  integer,public,parameter :: LOG_ALLREDUCE    = 40

  public :: timer_initialize
  public :: timer_set, timer_reset
  public :: timer_get, timer_thread_get
  public :: timer_begin, timer_end
  public :: timer_thread_begin, timer_thread_end
  public :: timer_show_hour, timer_show_min
  public :: timer_show_current_hour, timer_show_current_min

  public :: timer_reentrance_read, timer_reentrance_write
  public :: timer_write, timer_thread_write


  integer,private,parameter   :: LOG_SIZE = 50
  real(8),private,allocatable :: log_time(:)
  real(8),private,allocatable :: log_temp(:)

  real(8),private,allocatable :: log_time_t(:,:)
  real(8),private,allocatable :: log_temp_t(:,:)

private
contains
  subroutine timer_initialize
    use omp_lib, only: omp_get_max_threads
    implicit none
    allocate(log_time(0:LOG_SIZE - 1))
    allocate(log_temp(0:LOG_SIZE - 1))
    allocate(log_time_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    allocate(log_temp_t(0:LOG_SIZE - 1, 0:omp_get_max_threads()-1))
    call timer_reset
  end subroutine

  subroutine timer_set(e,t)
    implicit none
    integer,intent(in) :: e
    real(8),intent(in) :: t
    log_time(e) = t
    log_temp(e) = 0.d0
  end subroutine

  subroutine timer_reset(e)
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

  subroutine timer_reentrance_read(fd)
    implicit none
    integer,intent(in) :: fd
    read(fd) log_time(0:LOG_SIZE - 1)
    read(fd) log_temp(0:LOG_SIZE - 1)
  end subroutine

  subroutine timer_reentrance_write(fd)
    implicit none
    integer,intent(in) :: fd
    write(fd) log_time(0:LOG_SIZE - 1)
    write(fd) log_temp(0:LOG_SIZE - 1)
  end subroutine

  subroutine timer_begin(id)
    implicit none
    integer,intent(in) :: id
    log_temp(id) = get_wtime()
  end subroutine

  subroutine timer_end(id)
    implicit none
    integer,intent(in) :: id
    log_time(id) = log_time(id) + get_wtime() - log_temp(id)
  end subroutine

  subroutine timer_thread_begin(id)
    use omp_lib
    implicit none
    integer,intent(in) :: id
    integer :: tid
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timer_begin(id)
    end if
    log_temp_t(id,tid) = get_wtime()
  end subroutine

  subroutine timer_thread_end(id)
    use omp_lib
    implicit none
    integer,intent(in) :: id
    integer :: tid
    tid = omp_get_thread_num()
    if (tid == 0) then
      call timer_end(id)
    end if
    log_time_t(id,tid) = log_time_t(id,tid) + get_wtime() - log_temp_t(id,tid)
  end subroutine

  subroutine timer_show_hour(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,hour
    time = log_time(id)
    hour = time / 3600
    write(*,*) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timer_show_current_hour(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,hour
    time = get_wtime() - log_temp(id) + log_time(id)
    hour = time / 3600
    write(*,*) str,time,'sec =',hour,'hour'
  end subroutine

  subroutine timer_show_min(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,mini
    time = log_time(id)
    mini = time / 60
    write(*,*) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timer_show_current_min(str, id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: id
    real(8) :: time,mini
    time = get_wtime() - log_temp(id) + log_time(id)
    mini = time / 60
    write(*,*) str,time,'sec =',mini,'min'
  end subroutine

  subroutine timer_write(fd,str,id)
    implicit none
    character(*),intent(in) :: str
    integer,intent(in)      :: fd,id
    real(8) :: time
    time = log_time(id)
    write(fd,*) str,time,'sec'
  end subroutine

  subroutine timer_thread_write(fd,str,id)
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

  function timer_get(id)
    implicit none
    integer,intent(in) :: id
    real(8)            :: timer_get
    timer_get = log_time(id)
  end function

  function timer_thread_get(id,tid)
    implicit none
    integer,intent(in) :: id,tid
    real(8)            :: timer_thread_get
    timer_thread_get = log_time_t(id,tid)
  end function
end module
