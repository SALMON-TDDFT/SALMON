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
module misc_routines
  implicit none

  public :: floor_pow2, ceiling_pow2
  public :: gen_logfilename
  public :: get_wtime
  public :: create_directory

private
contains
  function floor_pow2(n)
    implicit none
    integer            :: floor_pow2
    integer,intent(in) :: n
    integer :: k

    k = 1
    do while(k < n)
      k = k * 2
    end do
    if (k /= n) k = ishft(k, -1)
    floor_pow2 = k
  end function

  function ceiling_pow2(n)
    implicit none
    integer            :: ceiling_pow2
    integer,intent(in) :: n
    integer :: k

    k = 1
    do while(k < n)
      k = k * 2
    end do
    ceiling_pow2 = k
  end function ceiling_pow2

  ! input  : base filename
  ! output : <base filename>_YYYYMMDD_hhmmss.log
  function gen_logfilename(filename)
    implicit none
    character(100)          :: gen_logfilename
    character(*),intent(in) :: filename
    character(8)  :: d
    character(10) :: t
    call date_and_time(date=d,time=t)

    write (gen_logfilename,'(A)') filename//'_'//d//'_'//t(1:6)//'.log'
  end function

  function get_wtime()
    implicit none
    real(8) :: get_wtime
    real(8) :: omp_get_wtime
    get_wtime = omp_get_wtime()
  end function

  ! NOTE: execute_command_line() is standardized at Fortran2008.
  !       In specification, `Execute command line` is defined this feature.
  subroutine create_directory(dirpath)
    implicit none
    character(*), intent(in) :: dirpath
    call execute_command_line('mkdir -p '//adjustl(trim(dirpath)), wait=.true.)
  end subroutine
end module
