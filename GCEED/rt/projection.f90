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
subroutine projection(tzpsi)
use salmon_parallel, only: nproc_group_grid, nproc_group_global, nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
use scf_data
use new_world_sub
use allocate_mat_sub
implicit none
complex(8) :: tzpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                    mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
integer :: ix,iy,iz,iob
integer :: iob_myob,icorr_p,job,job_allob
complex(8) :: coef_mat(itotMST,itotMST0,1,1)
complex(8) :: coef_mat2(itotMST,itotMST0,1,1)
real(8) :: coef(itotMST0,1,1)
complex(8) :: cbox
integer :: iobmax
integer :: iroot
complex(8),parameter :: zi=(0.d0,1.d0)
character(100) :: projOutFile

call calc_pmax(iobmax)


if(iSCFRT==2)then
  if(iwrite_projnum==1)then
    write(fileNumber, '(i8)') itt
    projOutFile = trim("proj.")//adjustl(fileNumber)
    open(61,file=projOutFile)
  end if
end if

coef_mat=0.d0

do iob=1,itotMST0
  call calc_myob(iob,iob_myob)
  call check_corrkob(iob,icorr_p)
  if(icorr_p==1)then
!$OMP parallel do private(iz,iy,ix) collapse(3)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      cmatbox_m(ix,iy,iz)=zpsi_t0(ix,iy,iz,iob_myob,1)
    end do
    end do
    end do
  end if
  call calc_iroot(iob,iroot)
  call comm_bcast(cmatbox_m,nproc_group_grid,iroot)
  do job=1,iobmax
    cbox=0.d0
!$OMP parallel do reduction(+:cbox) private(iz,iy,ix) collapse(3)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      cbox=cbox+conjg(tzpsi(ix,iy,iz,job,1))*cmatbox_m(ix,iy,iz)
    end do
    end do
    end do
    call calc_allob(job,job_allob)
    coef_mat(job_allob,iob,1,1)=coef_mat(job_allob,iob,1,1)+cbox
  end do
end do

call comm_summation(coef_mat,coef_mat2,itotMST*itotMST0,nproc_group_global)

coef=0.d0
do iob=1,itotMST0
  do job=1,itotMST
    coef(iob,1,1)=coef(iob,1,1)+abs(coef_mat2(job,iob,1,1)*Hvol)**2
  end do
end do
if(comm_is_root(nproc_id_global))then
  write(41,'(200f14.8)') dble(itt)*dt*2.41888d-2, &
  & (coef(iwrite_projection_ob(iob),iwrite_projection_k(iob),1),iob=1,num_projection),  &
    sum(coef(1:itotMST,:,1)),sum(coef(1:itotMST0,:,1))
end if
if(mod(itt,100)==0)then
  if(comm_is_root(nproc_id_global))then
    do iob=1,itotMST0
      write(*,'(a12,2i6,f16.8)') "projection",iob,coef(iob,1,1)
    end do
  end if
end if

if(iSCFRT==2)then
  if(iwrite_projnum==1)then
    close(61)
  end if
end if

end subroutine projection
