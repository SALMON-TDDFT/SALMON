!
!  Copyright 2018 SALMON developers
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
subroutine inner_product7(matbox1,matbox2,rbox2)
  use salmon_parallel, only: nproc_group_korbital
  use salmon_communication, only: comm_summation
  use misc_routines, only: get_wtime
  use scf_data
  use new_world_sub
  implicit none
  integer :: ix,iy,iz,iob,iob_allob
  real(8) :: matbox1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum)
  real(8) :: matbox2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum)
  real(8) :: rbox,rbox1(itotMST),rbox2(itotMST)
  
  rbox1(:)=0.d0
 
  do iob=1,iobnum
    call calc_allob(iob,iob_allob)
    rbox=0.d0
    !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : rbox)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      rbox=rbox+matbox1(ix,iy,iz,iob)*matbox2(ix,iy,iz,iob)
    end do
    end do
    end do
    rbox1(iob_allob)=rbox*Hvol
  end do
  
  elp3(186)=get_wtime()
  call comm_summation(rbox1,rbox2,itotMST,nproc_group_korbital)
  elp3(187)=get_wtime()
  elp3(190)=elp3(190)+elp3(187)-elp3(186)
  
end subroutine inner_product7
