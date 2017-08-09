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
subroutine read_input_gceed(cfunction2)
  use inputoutput
  use salmon_communication
  use salmon_parallel
  implicit none
  character(30),intent(out) :: cfunction2
  namelist / group_function2 / cfunction2

!  if(procid==0)then
!    open(fh_namelist, file='.namelist.tmp', status='old')
!    read(fh_namelist,nml=group_function2)
!    close(fh_namelist)
!  end if
!  call comm_bcast(cfunction2,nproc_group_global)

  cfunction2 = trim(calc_mode)

end subroutine read_input_gceed
