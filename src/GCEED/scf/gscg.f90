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
!======================================= Conjugate-Gradient minimization

subroutine sgscg(psi_in,iflag)
use salmon_parallel, only: nproc_group_grid, nproc_group_global, nproc_group_korbital
use salmon_communication, only: comm_summation, comm_bcast
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
use inner_product_sub
use hpsi2_sub
!$ use omp_lib
implicit none

real(8):: psi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),  &
               1:iobnum,1)
integer :: iter,iob,job,iflag
integer :: ix,iy,iz
integer :: is,iobsta(2),iobend(2)
real(8) :: sum0,sum1
real(8) :: sum_ob0(itotMST)
real(8) :: sum_obmat0(itotMST,itotMST),sum_obmat1(itotMST,itotMST)
real(8) :: xkHxk_ob(itotMST),xkxk_ob(itotMST),Rk_ob(itotMST)
real(8) :: gkgk_ob(itotMST),pkpk_ob(itotMST),xkpk_ob(itotMST)
real(8) :: pkHxk_ob(itotMST),pkHpk_ob(itotMST)
real(8) :: uk,alpha,Ak,Bk,Ck
real(8) , allocatable :: gk(:,:,:)
real(8) :: elp2(2000)
real(8):: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd)
real(8):: htpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
integer :: iob_myob,job_myob
integer :: iob_allob
integer :: icorr,jcorr,icorr_iob,icorr_job
integer :: iroot
integer :: is_sta,is_end

allocate (gk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

iwk_size=2
call make_iwksta_iwkend

call set_isstaend(is_sta,is_end)

!$OMP parallel do private(iz,iy,ix) collapse(2)
do iz=mg_sta(3)-Nd,mg_end(3)+Nd
do iy=mg_sta(2)-Nd,mg_end(2)+Nd
do ix=mg_sta(1)-Nd,mg_end(1)+Nd
  tpsi(ix,iy,iz)=0.d0
end do
end do
end do

elp2(:)=0d0
elp2(1)=get_wtime()

if(ilsda == 0)then
  iobsta(1)=1
  iobend(1)=itotMST
else if(ilsda == 1)then
  iobsta(1)=1
  iobend(1)=MST(1)
  iobsta(2)=MST(1)+1
  iobend(2)=itotMST
end if

do iob=1,iobnum
  call calc_allob(iob,iob_allob)

!$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rxk_ob(ix,iy,iz,iob)=psi_in(ix,iy,iz,iob,1)
    tpsi(ix,iy,iz)=rxk_ob(ix,iy,iz,iob)
  end do
  end do
  end do
  call hpsi2(tpsi,htpsi,iob_allob,1,0,0)
  
!$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhxk_ob(ix,iy,iz,iob)=htpsi(ix,iy,iz)
  end do
  end do
  end do
end do

call inner_product7(rxk_ob,rhxk_ob,xkHxk_ob)

xkxk_ob(:)=1.d0 
Rk_ob(:)=xkHxk_ob(:)/xkxk_ob(:)

Iteration : do iter=1,Ncg
elp2(2)=get_wtime()
  do iob=1,iobnum
    call calc_allob(iob,iob_allob)
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      rgk_ob(ix,iy,iz,iob) = 2*( rhxk_ob(ix,iy,iz,iob) - Rk_ob(iob_allob)*rxk_ob(ix,iy,iz,iob) )
    end do
    end do
    end do
  end do
 if(nproc_ob==1)then
    do is=is_sta,is_end
    do iob=iobsta(is),iobend(is)
      do job=iobsta(is),iob-1
        sum0=0.d0
  !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          sum0=sum0+psi_in(ix,iy,iz,job,1)*rgk_ob(ix,iy,iz,iob)
        end do
        end do
        end do
        sum_obmat0(iob,job)=sum0*Hvol
      end do
    end do 
    end do
    call comm_summation(sum_obmat0,sum_obmat1,itotMST*itotMST,nproc_group_global)
    do is=is_sta,is_end
    do iob=iobsta(is),iobend(is)
      do job=iobsta(is),iob-1
  !$omp parallel do private(iz,iy,ix) collapse(2)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rgk_ob(ix,iy,iz,iob)=rgk_ob(ix,iy,iz,iob)-sum_obmat1(iob,job)*psi_in(ix,iy,iz,job,1)
        end do
        end do
        end do
      end do
    end do
    end do
  else
    do iob=iobsta(is),iobend(is)
      call calc_myob(iob,iob_myob)
      call check_corrkob(iob,1,icorr_iob)
      do job=iobsta(is),iob-1
        call calc_myob(job,job_myob)
        call check_corrkob(job,1,icorr_job)
        if(icorr_job==1)then
  !$omp parallel do private(iz,iy,ix) collapse(2)
          do iz=mg_sta(3),mg_end(3)
          do iy=mg_sta(2),mg_end(2)
          do ix=mg_sta(1),mg_end(1)
            matbox_m(ix,iy,iz)=psi_in(ix,iy,iz,job_myob,1)
          end do
          end do
          end do
        end if
        call calc_iroot(job,iroot)
        call comm_bcast(matbox_m,nproc_group_grid,iroot)
        sum0=0.d0
  !$omp parallel do private(iz,iy,ix) collapse(2) reduction(+ : sum0)
        do iz=iwk3sta(3),iwk3end(3)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          sum0=sum0+matbox_m(ix,iy,iz)*rgk_ob(ix,iy,iz,iob_myob)
        end do
        end do
        end do
        sum0=sum0*Hvol
        call comm_summation(sum0,sum1,nproc_group_korbital)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          rgk_ob(ix,iy,iz,iob_myob)=rgk_ob(ix,iy,iz,iob_myob)-sum1*matbox_m(ix,iy,iz)
        end do
        end do
        end do
      end do
    end do
  end if 
 call inner_product7(rgk_ob,rgk_ob,sum_ob0)
 if ( iter==1 ) then
    do iob=1,iobnum
      call calc_allob(iob,iob_allob)
!$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rpk_ob(ix,iy,iz,iob) = -rgk_ob(ix,iy,iz,iob)
      end do
      end do
      end do
    end do
  else
    do iob=1,iobnum
      call calc_allob(iob,iob_allob)
      uk=sum_ob0(iob_allob)/gkgk_ob(iob_allob)
!$OMP parallel do private(iz,iy,ix)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rpk_ob(ix,iy,iz,iob) = -rgk_ob(ix,iy,iz,iob) + uk*rpk_ob(ix,iy,iz,iob)
      end do
      end do
      end do
    end do
  end if 
  gkgk_ob(:)=sum_ob0(:)
  call inner_product7(rxk_ob,rpk_ob,xkpk_ob)
  call inner_product7(rpk_ob,rpk_ob,pkpk_ob)
  call inner_product7(rpk_ob,rhxk_ob,pkHxk_ob)

  do iob=1,iobnum
    call calc_allob(iob,iob_allob)
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      tpsi(ix,iy,iz) = rpk_ob(ix,iy,iz,iob)
    end do
    end do
    end do
    call hpsi2(tpsi,htpsi,iob_allob,1,0,0)
    
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
       rgk_ob(ix,iy,iz,iob)=htpsi(ix,iy,iz)
    end do
    end do
    end do
  end do
  call inner_product7(rpk_ob,rgk_ob,pkHpk_ob)
 do iob=1,iobnum
    call calc_allob(iob,iob_allob)
    Ak=pkHpk_ob(iob_allob)*xkpk_ob(iob_allob)-pkHxk_ob(iob_allob)*pkpk_ob(iob_allob)
    Bk=pkHpk_ob(iob_allob)*xkxk_ob(iob_allob)-xkHxk_ob(iob_allob)*pkpk_ob(iob_allob)
    Ck=pkHxk_ob(iob_allob)*xkxk_ob(iob_allob)-xkHxk_ob(iob_allob)*xkpk_ob(iob_allob)
    alpha=(-Bk+sqrt(Bk*Bk-4*Ak*Ck))/(2*Ak)

!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      rxk_ob(ix,iy,iz,iob) = rxk_ob(ix,iy,iz,iob) + alpha*rpk_ob(ix,iy,iz,iob)
      rhxk_ob(ix,iy,iz,iob) = rhxk_ob(ix,iy,iz,iob) + alpha*rgk_ob(ix,iy,iz,iob)
    end do
    end do
    end do
  end do
  call inner_product7(rxk_ob,rhxk_ob,xkHxk_ob)
  call inner_product7(rxk_ob,rxk_ob,xkxk_ob)
  Rk_ob(:)=xkHxk_ob(:)/xkxk_ob(:)


end do Iteration

call inner_product7(rxk_ob,rxk_ob,sum_ob0)
do iob=1,iobnum
  call calc_allob(iob,iob_allob)
!$OMP parallel do private(iz,iy,ix) collapse(2)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    psi_in(ix,iy,iz,iob,1)=rxk_ob(ix,iy,iz,iob)/sqrt(sum_ob0(iob_allob))
  end do
  end do
  end do
end do

if(iflag.eq.1) then
  iflag=0
end if

deallocate (gk)
return

end subroutine sgscg
