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
subroutine exc_cor_ns
  use salmon_parallel, only: nproc_group_h
  use salmon_communication, only: comm_summation
  use salmon_xc, only: calc_xc
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  use sendrecvh_sub
  implicit none
  integer :: ix,iy,iz,is
  real(8) :: tot_exc
  real(8),allocatable :: rhd(:,:,:), delr(:,:,:,:)
  integer :: iwk_dum

  if(ilsda==0)then
    do iz=1,ng_num(3)
    do iy=1,ng_num(2)
    do ix=1,ng_num(1)
      rho_tmp(ix,iy,iz)=rho(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1)
    end do
    end do
    end do
  else if(ilsda==1)then
    do is=1,2
    do iz=1,ng_num(3)
    do iy=1,ng_num(2)
    do ix=1,ng_num(1)
      rho_s_tmp(ix,iy,iz,is)=rho_s(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1,is)
    end do
    end do
    end do
    end do
  end if

  if(xc=='pz'.or.xc=='PZ')then
    continue
  else
    allocate (rhd (ng_sta(1)-Ndh:ng_end(1)+Ndh, &
                   ng_sta(2)-Ndh:ng_end(2)+Ndh, &
                   ng_sta(3)-Ndh:ng_end(3)+Ndh))
    allocate (delr(ng_sta(1):ng_end(1), &
                   ng_sta(2):ng_end(2), &
                   ng_sta(3):ng_end(3),3))
  
!$OMP parallel do private(ix,iy,iz)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      rhd(ix,iy,iz)=dble(rho(ix,iy,iz))
    enddo
    enddo
    enddo
  
!$omp end parallel do
    iwk_dum=iwk_size
    iwk_size=12
    call make_iwksta_iwkend
    call sendrecvh(rhd)
    iwk_size=iwk_dum

!$OMP parallel do private(ix,iy,iz)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
       delr(ix,iy,iz,1)= (1.0d0/Hgs(1))*(   &
          bN1*( rhd(ix+1,iy,iz) - rhd(ix-1,iy,iz))  &
        + bN2*( rhd(ix+2,iy,iz) - rhd(ix-2,iy,iz)) &
        + bN3*( rhd(ix+3,iy,iz) - rhd(ix-3,iy,iz)) &
        + bN4*( rhd(ix+4,iy,iz) - rhd(ix-4,iy,iz)))
       delr(ix,iy,iz,2)= (1.0d0/Hgs(2))*( &
          bN1*( rhd(ix,iy+1,iz) - rhd(ix,iy-1,iz)) &
        + bN2*( rhd(ix,iy+2,iz) - rhd(ix,iy-2,iz)) &
        + bN3*( rhd(ix,iy+3,iz) - rhd(ix,iy-3,iz)) &
        + bN4*( rhd(ix,iy+4,iz) - rhd(ix,iy-4,iz)))
       delr(ix,iy,iz,3)= (1.0d0/Hgs(3))*( &
          bN1*( rhd(ix,iy,iz+1) - rhd(ix,iy,iz-1)) &
        + bN2*( rhd(ix,iy,iz+2) - rhd(ix,iy,iz-2)) &
        + bN3*( rhd(ix,iy,iz+3) - rhd(ix,iy,iz-3)) &
        + bN4*( rhd(ix,iy,iz+4) - rhd(ix,iy,iz-4)))
    enddo
    enddo
    enddo
!$omp end parallel do
  end if

  if(xc=='pz'.or.xc=='PZ')then
    if(ilsda==0)then
      call calc_xc(xc_func, rho=rho_tmp, eexc=eexc_tmp, vxc=vxc_tmp)
    else if(ilsda==1)then
      call calc_xc(xc_func, rho_s=rho_s_tmp, eexc=eexc_tmp, vxc_s=vxc_s_tmp)
    end if
  else
    call calc_xc(xc_func, rho=rho_tmp, grho=delr, eexc=eexc_tmp, vxc=vxc_tmp)
  end if

  if(ilsda==0)then
    do iz=1,ng_num(3)
    do iy=1,ng_num(2)
    do ix=1,ng_num(1)
      vxc(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1)=vxc_tmp(ix,iy,iz)
    end do
    end do
    end do
  else if(ilsda==1)then 
    do is=1,2
    do iz=1,ng_num(3)
    do iy=1,ng_num(2)
    do ix=1,ng_num(1)
      vxc_s(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1,is)=vxc_s_tmp(ix,iy,iz,is)
    end do
    end do
    end do
    end do
  end if

  tot_exc=0.d0
  do iz=1,ng_num(3)
  do iy=1,ng_num(2)
  do ix=1,ng_num(1)
    tot_exc=tot_exc+eexc_tmp(ix,iy,iz)
  end do
  end do
  end do
  tot_exc=tot_exc*hvol
 
  call comm_summation(tot_exc,Exc,nproc_group_h)

  return
end subroutine exc_cor_ns
