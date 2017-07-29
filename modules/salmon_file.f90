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
!-----------------------------------------------------------------------------------------
module salmon_file
  implicit none
  private

  integer, parameter :: fh_start = 1000

  public :: get_filehandle
contains
!--------------------------------------------------------------------------------
!! Return a unit number available to open file
  integer function get_filehandle() result(fh)
    implicit none
    logical :: flag
    
    fh = fh_start - 1
    flag = .true.
    do while (flag)
      fh = fh + 1
      inquire(unit=fh, opened=flag)
    end do
    return
  end function get_filehandle
  
end module salmon_file
!--------------------------------------------------------------------------------
