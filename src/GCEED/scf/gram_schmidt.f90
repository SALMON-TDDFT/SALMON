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
!======================================== Gram-Schmidt orthogonalization
SUBROUTINE Gram_Schmidt_ns
use salmon_parallel, only: nproc_group_kgrid, nproc_group_global
use salmon_communication, only: comm_summation, comm_bcast
use scf_data
use new_world_sub
use allocate_mat_sub
implicit none

integer :: iob,job,iob_myob,job_allob
integer :: is,pstart(2),pend(2)
integer :: ix,iy,iz
real(8) :: ovrp(1:itotMST),ovrp2(1:itotMST)
real(8) :: rbox,rbox2
integer :: iroot
integer :: icorr_p
integer :: is_sta,is_end

call set_isstaend(is_sta,is_end)

if(ilsda == 0)then
  pstart(1)=1
  pend(1)=itotMST
else if(ilsda == 1)then
  pstart(1)=1
  pend(1)=MST(1)
  pstart(2)=MST(1)+1
  pend(2)=itotMST
end if

do is=is_sta,is_end
do iob=pstart(is),pend(is)
  call calc_myob(iob,iob_myob)
  call check_corrkob(iob,1,icorr_p)
  if(icorr_p==1)then
!$OMP parallel do private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      matbox_m(ix,iy,iz)=psi(ix,iy,iz,iob_myob,1)
    end do
    end do
    end do
  end if
  call calc_iroot(iob,iroot)
  call comm_bcast(matbox_m,nproc_group_kgrid,iroot)

  ovrp=0.d0
  do job=1,iobnum
    call calc_allob(job,job_allob)
    if(job_allob >= pstart(is) .and. job_allob <= iob-1)then
      rbox=0.d0
!$OMP parallel do reduction ( + : rbox ) private(iz,iy,ix)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rbox=rbox+psi(ix,iy,iz,job,1)*matbox_m(ix,iy,iz)*Hvol
      end do
      end do
      end do
      ovrp(job_allob)=ovrp(job_allob)+rbox
    end if
  end do

  call comm_summation(ovrp,ovrp2,itotMST,nproc_group_global)

  matbox_m=0.d0
  do job=1,iobnum
    call calc_allob(job,job_allob)
    if(job_allob >= pstart(is) .and. job_allob <= iob-1)then
!$OMP parallel do private(iz,iy,ix)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        matbox_m(ix,iy,iz)=matbox_m(ix,iy,iz)-ovrp2(job_allob)*psi(ix,iy,iz,job,1)
      end do
      end do
      end do
    end if
  end do

  call comm_summation(matbox_m,matbox_m2,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_kgrid)

  rbox=0.d0
  call check_corrkob(iob,1,icorr_p)
  if(icorr_p==1)then
!$OMP parallel do reduction ( + : rbox ) private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      psi(ix,iy,iz,iob_myob,1)=psi(ix,iy,iz,iob_myob,1)+matbox_m2(ix,iy,iz)
      rbox=rbox+psi(ix,iy,iz,iob_myob,1)**2
    end do
    end do
    end do
  end if

  call comm_summation(rbox,rbox2,nproc_group_global)

  if(icorr_p==1)then
!$OMP parallel do private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      psi(ix,iy,iz,iob_myob,1)=psi(ix,iy,iz,iob_myob,1)/sqrt(rbox2*Hvol)
    end do
    end do
    end do
  end if

end do
end do

return

END SUBROUTINE Gram_Schmidt_ns
