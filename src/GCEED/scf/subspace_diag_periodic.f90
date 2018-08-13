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
subroutine subspace_diag_periodic

use salmon_parallel, only: nproc_group_korbital, nproc_group_k, nproc_group_kgrid
use salmon_communication, only: comm_bcast, comm_summation
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
use hpsi2_sub
implicit none
integer :: ii,jj,ik
integer :: ix,iy,iz
complex(8),allocatable :: Amat(:,:)
complex(8),allocatable :: Amat2(:,:)
complex(8),allocatable :: Smat(:,:)
complex(8):: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd)
complex(8),allocatable :: htpsi(:,:,:)
complex(8),allocatable :: htpsi_groupob(:,:,:,:)
complex(8),allocatable :: ztpsi_groupob(:,:,:,:)
complex(8) :: cbox
real(8) :: rbox,rbox1
complex(8),allocatable :: evec(:,:)
integer :: iter,ier2
integer :: iroot
integer :: j_myob,i_allob,j_allob
integer :: icorr_j
integer :: is,is_sta,is_end
integer :: iobsta(2),iobend(2)

elp3(301)=get_wtime()

allocate(htpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(htpsi_groupob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
allocate(ztpsi_groupob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))

iwk_size=2
call make_iwksta_iwkend

if(ilsda == 0)then
  iobsta(1)=1
  iobend(1)=itotMST
else if(ilsda == 1)then
  iobsta(1)=1
  iobend(1)=MST(1)
  iobsta(2)=MST(1)+1
  iobend(2)=itotMST
end if

!$OMP parallel do private(iz,iy,ix)
do iz=mg_sta(3)-Nd,mg_end(3)+Nd
do iy=mg_sta(2)-Nd,mg_end(2)+Nd
do ix=mg_sta(1)-Nd,mg_end(1)+Nd
  tpsi(ix,iy,iz)=0.d0
end do
end do
end do

elp3(302)=get_wtime()
elp3(352)=elp3(352)+elp3(302)-elp3(301)

elp3(303)=get_wtime()
elp3(353)=elp3(353)+elp3(303)-elp3(302)

call set_isstaend(is_sta,is_end)

do ik=k_sta,k_end
do is=is_sta,is_end

  if(ifMST(is)>=1.and.MST(is)>=1)then

    iter=iobend(is)-iobsta(is)+1
  
    allocate(evec(iter,iter))
    allocate(Amat(iter,iter))
    allocate(Amat2(iter,iter))
    allocate(Smat(iter,iter))
  
!$OMP parallel do private(jj,ii)
    do jj=1,iter
      do ii=1,iter
        Amat(ii,jj)=0.d0
      end do
    end do
  
!do jj=1,itotMST
    do jj=1,iobnum
      call calc_allob(jj,j_allob)
      if(j_allob>=iobsta(is).and.j_allob<=iobend(is))then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          tpsi(ix,iy,iz)=zpsi(ix,iy,iz,jj,ik)
        end do
        end do
        end do
        call hpsi2(tpsi,htpsi,j_allob,ik,0,0)
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          htpsi_groupob(ix,iy,iz,jj)=htpsi(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do
  
    do jj=iobsta(is),iobend(is)
      call calc_myob(jj,j_myob)
      call check_corrkob(jj,ik,icorr_j)
      if(icorr_j==1)then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          htpsi(ix,iy,iz)=htpsi_groupob(ix,iy,iz,j_myob)
        end do
        end do
        end do
      end if
      call calc_iroot(jj,iroot)
      call comm_bcast(htpsi,nproc_group_kgrid,iroot)
      do ii=1,iobnum
        call calc_allob(ii,i_allob)
        if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
          cbox=0.d0
!$OMP parallel do private(iz,iy,ix) reduction (+ : cbox)
          do iz=mg_sta(3),mg_end(3)
          do iy=mg_sta(2),mg_end(2)
          do ix=mg_sta(1),mg_end(1)
             cbox=cbox+conjg(zpsi(ix,iy,iz,ii,ik))*htpsi(ix,iy,iz)
          end do
          end do
          end do
          Amat(i_allob-iobsta(is)+1,jj-iobsta(is)+1)=cbox*Hvol
        end if
      end do
    end do
  
    call comm_summation(Amat,Amat2,iter*iter,nproc_group_k)
  
    call eigen_subdiag_periodic(Amat2,evec,iter,ier2)
  
    do jj=1,iobnum
      call calc_allob(jj,j_allob)
      if(j_allob>=iobsta(is).and.j_allob<=iobend(is))then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          ztpsi_groupob(ix,iy,iz,jj)=zpsi(ix,iy,iz,jj,ik)
          zpsi(ix,iy,iz,jj,ik)=0.d0
        end do
        end do
        end do
      end if
    end do
    
    do jj=iobsta(is),iobend(is)
      call calc_myob(jj,j_myob)
      call check_corrkob(jj,ik,icorr_j)
      if(icorr_j==1)then
!$OMP parallel do private(iz,iy,ix)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          htpsi(ix,iy,iz)=ztpsi_groupob(ix,iy,iz,j_myob) ! htpsi is making a role of original tpsi
        end do
        end do
        end do
      end if
      call calc_iroot(jj,iroot)
      call comm_bcast(htpsi,nproc_group_kgrid,iroot)
      do ii=1,iobnum
        call calc_allob(ii,i_allob)
        if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
!$OMP parallel do private(iz,iy,ix)
          do iz=mg_sta(3),mg_end(3)
          do iy=mg_sta(2),mg_end(2)
          do ix=mg_sta(1),mg_end(1)
            zpsi(ix,iy,iz,ii,ik)=zpsi(ix,iy,iz,ii,ik)+evec(jj-iobsta(is)+1,i_allob-iobsta(is)+1)*htpsi(ix,iy,iz)
          end do
          end do
          end do
        end if
      end do
    end do
  
    do ii=1,iobnum
      call calc_allob(ii,i_allob)
      if(i_allob>=iobsta(is).and.i_allob<=iobend(is))then
        rbox=0.d0
!$OMP parallel do private(iz,iy,ix) reduction(+:rbox)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rbox=rbox+abs(zpsi(ix,iy,iz,ii,ik))**2
        end do
        end do
        end do
        call comm_summation(rbox,rbox1,nproc_group_korbital)
!$OMP parallel do private(iz,iy,ix) 
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          zpsi(ix,iy,iz,ii,ik)=zpsi(ix,iy,iz,ii,ik)/sqrt(rbox1*Hvol)
        end do
        end do
        end do
      end if
    end do
  
    deallocate(Amat,Amat2,Smat)
    deallocate(evec)

  end if

end do
end do

deallocate(htpsi,htpsi_groupob,ztpsi_groupob)

end subroutine subspace_diag_periodic
