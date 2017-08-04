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
subroutine writeestatic
  use salmon_parallel, only: nproc_group_global
  use salmon_communication, only: comm_summation
  use scf_data
  use allocate_mat_sub
  implicit none
  integer :: ix,iy,iz,jj
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum
  character(20) :: header_unit

  do jj=1,3

    !$OMP parallel do collapse(2) private(iz,iy,ix)
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do
    
    if(jj==1)then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        matbox_l(ix,iy,iz)=Ex_static(ix,iy,iz)
      end do
      end do
      end do
    else if(jj==2)then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        matbox_l(ix,iy,iz)=Ey_static(ix,iy,iz)
      end do
      end do
      end do
    else if(jj==3)then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        matbox_l(ix,iy,iz)=Ez_static(ix,iy,iz)
      end do
      end do
      end do
    end if
   
    if(format3d=='avs')then
      !$OMP parallel do collapse(2) private(iz,iy,ix)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        matbox_l(ix,iy,iz)=matbox_l(ix,iy,iz)*5.14223d1
      end do
      end do
      end do
    end if
 
    call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
  
    write(filenum, '(i8)') itt
    if(jj==1)then
      suffix = "Exsta"//adjustl(filenum)
      phys_quantity = "exsta"
    else if(jj==2)then
      suffix = "Eysta"//adjustl(filenum)
      phys_quantity = "eysta"
    else if(jj==3)then
      suffix = "Ezsta"//adjustl(filenum)
      phys_quantity = "ezsta"
    end if

    if(format3d=='avs')then
      header_unit = "V/A"
      call writeavs(103,suffix,header_unit,matbox_l2)
    else if(format3d=='cube')then
      call writecube(103,suffix,phys_quantity,matbox_l2)
    end if

  end do  
 
end subroutine writeestatic
