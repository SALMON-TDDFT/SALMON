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
!======================================== Gram-Schmidt orthogonalization
subroutine Gram_Schmidt_periodic
use salmon_parallel, only: nproc_group_kgrid, nproc_group_korbital, nproc_group_k
use salmon_communication, only: comm_summation, comm_bcast
use scf_data
use new_world_sub
implicit none

integer :: iob,q,iob_myob,q_allob
integer :: is,iobsta(2),iobend(2)
integer :: ik
integer :: ix,iy,iz
real(8) :: rbox,rbox2
complex(8) :: zovrp(1:itotMST),zovrp2(1:itotMST)
complex(8) :: cbox
complex(8),allocatable :: cmatbox(:,:,:),cmatbox2(:,:,:)
integer :: icorr_p
integer :: iroot
integer :: is_sta,is_end

call set_isstaend(is_sta,is_end)

if(ilsda == 0)then
  iobsta(1)=1
  iobend(1)=itotMST
else if(ilsda == 1)then
  iobsta(1)=1
  iobend(1)=MST(1)
  iobsta(2)=MST(1)+1
  iobend(2)=itotMST
end if

allocate(cmatbox(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(cmatbox2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

do ik=k_sta,k_end
  do is=is_sta,is_end
  do iob=iobsta(is),iobend(is)
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,ik,icorr_p)
    if(icorr_p==1)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        cmatbox(ix,iy,iz)=zpsi(ix,iy,iz,iob_myob,ik)
      end do
      end do
      end do
    end if
    call calc_iroot(iob,iroot)
    call comm_bcast(cmatbox,nproc_group_kgrid,iroot)

    zovrp=0.d0
    do q=1,iobnum
      call calc_allob(q,q_allob)
      if(q_allob >= iobsta(is) .and. q_allob <= iob-1)then
        cbox=0.d0
!$OMP parallel do private(iz,iy,ix) collapse(2) reduction ( + : cbox )
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cbox=cbox+conjg(zpsi(ix,iy,iz,q,ik))*cmatbox(ix,iy,iz)*Hvol
        end do
        end do
        end do
        zovrp(q_allob)=zovrp(q_allob)+cbox
      end if
    end do

    call comm_summation(zovrp,zovrp2,itotMST,nproc_group_k)

    cmatbox=0.d0
    do q=1,iobnum
      call calc_allob(q,q_allob)
      if(q_allob >= iobsta(is) .and. q_allob <= iob-1)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox(ix,iy,iz)=cmatbox(ix,iy,iz)-zovrp2(q_allob)*zpsi(ix,iy,iz,q,ik)
        end do
        end do
        end do
      end if
    end do

    call comm_summation(cmatbox,cmatbox2,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_kgrid)

    rbox=0.d0
    if(icorr_p==1)then
!$OMP parallel do private(iz,iy,ix) collapse(2) reduction ( + : rbox )
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        zpsi(ix,iy,iz,iob_myob,ik)=zpsi(ix,iy,iz,iob_myob,ik)+cmatbox2(ix,iy,iz)
        rbox=rbox+abs(zpsi(ix,iy,iz,iob_myob,ik))**2
      end do
      end do
      end do
    end if

    call comm_summation(rbox,rbox2,nproc_group_korbital)

    if(icorr_p==1)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        zpsi(ix,iy,iz,iob_myob,ik)=zpsi(ix,iy,iz,iob_myob,ik)/sqrt(rbox2*Hvol)
      end do
      end do
      end do
    end if

  end do
  end do

end do

deallocate (cmatbox,cmatbox2)

return

end subroutine Gram_Schmidt_periodic
