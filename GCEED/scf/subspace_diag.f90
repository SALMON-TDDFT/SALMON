! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine subspace_diag

use scf_data
use inner_product_sub
use hpsi2_sub
use copy_psi_mesh_sub
implicit none
integer :: iob,job,ii,jj
integer :: jmax
integer :: ix,iy,iz,is
real(8),allocatable :: Amat(:,:)
real(8),allocatable :: Amat2(:,:)
real(8),allocatable :: Smat(:,:)
real(8),allocatable :: tpsi(:,:,:),htpsi(:,:,:)
real(8),allocatable :: psi_box(:,:,:,:)
real(8) :: rbox,rbox1
real(8),allocatable :: alpha(:),betav(:),eval(:)
real(8),allocatable :: evec(:,:),evec_box(:)
real(8) :: vmax
integer :: ier2
real(8) :: rbox_array(itotMST)
real(8) :: rbox_array2(itotMST)
integer :: is_sta,is_end
integer :: job_myob,iroot,icorr_j,iob_allob,job_allob
integer :: iter
integer :: iobsta(2),iobend(2)

elp3(301)=MPI_Wtime()

allocate(tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd))
allocate(htpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(psi_box(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))

iwk_size=2
call make_iwksta_iwkend

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

!$OMP parallel do
do iz=mg_sta(3)-Nd,mg_end(3)+Nd
do iy=mg_sta(2)-Nd,mg_end(2)+Nd
do ix=mg_sta(1)-Nd,mg_end(1)+Nd
  tpsi(ix,iy,iz)=0.d0
end do
end do
end do

do is=is_sta,is_end

  iter=iobend(is)-iobsta(is)+1
  allocate(evec(iter,iter))
  allocate(Amat(iter,iter))
  allocate(Amat2(iter,iter))
  allocate(Smat(iter,iter))

  do jj=1,iter
!$OMP parallel do
    do ii=1,iter
      Amat2(ii,jj)=0.d0
    end do
  end do

  do job=iobsta(is),iobend(is)
    call calc_myob(job,job_myob)
    call check_corrkob(job,icorr_j)
    if(icorr_j==1)then
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        tpsi(ix,iy,iz)=psi(ix,iy,iz,job_myob,1)
      end do
      end do
      end do
      call hpsi2(tpsi,htpsi,job,0,0)
    end if
    call calc_iroot(job,iroot)
    call MPI_Bcast(htpsi,mg_num(1)*mg_num(2)*mg_num(3),MPI_DOUBLE_PRECISION,iroot,newworld_comm_grid,ierr)
    
    do iob=1,iobnum
      call calc_allob(iob,iob_allob)
      if(iob_allob>=iobsta(is).and.iob_allob<=iobend(is))then
        rbox=0.d0
!$OMP parallel do reduction(+:rbox)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rbox=rbox+psi(ix,iy,iz,iob,1)*htpsi(ix,iy,iz)
        end do
        end do
        end do
        Amat2(iob_allob-iobsta(is)+1,job-iobsta(is)+1)=rbox*Hvol
      end if
    end do
  end do
  
  call MPI_Allreduce(Amat2,Amat,iter*iter,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  call eigen_subdiag(Amat,evec,iter,ier2)
 
  do job=1,iobnum
    call calc_allob(job,job_allob)
    if(job_allob>=iobsta(is).and.job_allob<=iobend(is))then
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        psi_box(ix,iy,iz,job)=psi(ix,iy,iz,job,1)
        psi(ix,iy,iz,job,1)=0.d0
      end do
      end do
      end do
    end if
  end do
   
  do job=iobsta(is),iobend(is)
    call calc_myob(job,job_myob)
    call check_corrkob(job,icorr_j)
    if(icorr_j==1)then
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        matbox_m(ix,iy,iz)=psi_box(ix,iy,iz,job_myob)
      end do
      end do
      end do
    end if
    call calc_iroot(job,iroot)
    call MPI_Bcast(matbox_m,mg_num(1)*mg_num(2)*mg_num(3),MPI_DOUBLE_PRECISION,iroot,newworld_comm_grid,ierr)
    do iob=1,iobnum
      call calc_allob(iob,iob_allob)
      if(iob_allob>=iobsta(is).and.iob_allob<=iobend(is))then
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          psi(ix,iy,iz,iob,1)=psi(ix,iy,iz,iob,1)+evec(job-iobsta(is)+1,iob_allob-iobsta(is)+1)*matbox_m(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do
  end do

  do iob=1,iobnum
    call calc_allob(iob,iob_allob)
    if(iob_allob>=iobsta(is).and.iob_allob<=iobend(is))then
      rbox=0.d0
!$OMP parallel do reduction(+:rbox)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rbox=rbox+abs(psi(ix,iy,iz,iob,1))**2
      end do
      end do
      end do
      call MPI_Allreduce(rbox,rbox1,1,MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_orbital,ierr)
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        psi(ix,iy,iz,iob,1)=psi(ix,iy,iz,iob,1)/sqrt(rbox1*Hvol)
      end do
      end do
      end do
    end if
  end do
  deallocate(Amat,Amat2,Smat)
  deallocate(evec)

end do

deallocate(htpsi,psi_box)

end subroutine subspace_diag
