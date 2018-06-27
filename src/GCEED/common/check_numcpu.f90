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
subroutine check_numcpu
use salmon_parallel, only: nproc_size_global
use scf_data
use new_world_sub
implicit none
integer :: j

if(nproc_k*nproc_ob*nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin(3)/=nproc_size_global)then
  write(*,*) "inumcpu_check error!"
  write(*,*) "number of cpu is not correct!"
  stop
end if
do j=1,3
  if(nproc_Mxin_s(j)<nproc_Mxin(j))then
    write(*,*) "inumcpu_check error!"
    write(*,*) "nproc_Mxin_s is smaller than nproc_Mxin."
    stop
  end if
end do
if(nproc_Mxin_s(1)*nproc_Mxin_s(2)*nproc_Mxin_s(3)>nproc_size_global)then
  write(*,*) "inumcpu_check error!"
  write(*,*) "product of nproc_Mxin_s is larger than nproc."
  stop
end if
if(mod(nproc_Mxin_s(1),nproc_Mxin(1))/=0)then
  write(*,*) "inumcpu_check error!"
  write(*,*) "nproc_Mxin_s(1) is not mutiple of nproc_Mxin(1)."
  stop
end if
if(mod(nproc_Mxin_s(2),nproc_Mxin(2))/=0)then
  write(*,*) "inumcpu_check error!"
  write(*,*) "nproc_Mxin_s(2) is not mutiple of nproc_Mxin(2)."
  stop
end if
if(mod(nproc_Mxin_s(3),nproc_Mxin(3))/=0)then
  write(*,*) "inumcpu_check error!"
  write(*,*) "nproc_Mxin_s(3) is not mutiple of nproc_Mxin(3)."
  stop
end if
nproc_Mxin_s_dm(1:3)=nproc_Mxin_s(1:3)/nproc_Mxin(1:3)

if(ilsda==1)then
  if(nproc_ob>1)then
    write(*,*) "Orbital parallelization is not currently supported. It will be supported in future."
    stop
  end if
end if

end subroutine check_numcpu
