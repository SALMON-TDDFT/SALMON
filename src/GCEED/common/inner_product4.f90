!
!  Copyright 2017-2019 SALMON developers
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
subroutine inner_product4(matbox1,matbox2,cbox2)
use salmon_parallel, only: nproc_group_korbital
use salmon_communication, only: comm_summation
use scf_data
use new_world_sub
!$ use omp_lib
implicit none
integer :: ix,iy,iz
complex(8) :: matbox1(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
complex(8) :: matbox2(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
complex(8) :: cbox,cbox2

cbox=0.d0
!$omp parallel do private(iz,iy,ix) reduction(+ : cbox)
do iz=iwk3sta(3),iwk3end(3)
do iy=iwk3sta(2),iwk3end(2)
do ix=iwk3sta(1),iwk3end(1)
  cbox=cbox+conjg(matbox1(ix,iy,iz))*matbox2(ix,iy,iz)
end do
end do
end do
cbox=cbox*Hvol
call comm_summation(cbox,cbox2,nproc_group_korbital)

end subroutine inner_product4
