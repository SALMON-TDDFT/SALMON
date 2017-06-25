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
!=======================================================================
subroutine inner_product3(matbox1,matbox2,rbox2)
use salmon_parallel, only: nproc_group_orbital
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
!$ use omp_lib
implicit none
integer :: ix,iy,iz
real(8) :: matbox1(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
real(8) :: matbox2(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
real(8) :: rbox,rbox2

rbox=0.d0
!$omp parallel do reduction(+ : rbox) collapse(3) private(iz,iy,ix)
do iz=iwk3sta(3),iwk3end(3)
do iy=iwk3sta(2),iwk3end(2)
do ix=iwk3sta(1),iwk3end(1)
  rbox=rbox+matbox1(ix,iy,iz)*matbox2(ix,iy,iz)
end do
end do
end do

elp3(186)=get_wtime()
call comm_summation(rbox,rbox2,nproc_group_orbital)
elp3(187)=get_wtime()
elp3(190)=elp3(190)+elp3(187)-elp3(186)

end subroutine inner_product3
