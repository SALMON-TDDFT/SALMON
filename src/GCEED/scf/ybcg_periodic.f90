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
!======================================= Conjugate-Gradient minimization

SUBROUTINE DTcg_periodic(psi_in,iflag)
use salmon_parallel, only: nproc_group_kgrid, nproc_group_korbital
use salmon_communication, only: comm_bcast, comm_summation
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
use inner_product_sub
use allocate_mat_sub
use hpsi2_sub
!$ use omp_lib
implicit none

complex(8) :: psi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),  &
               1:iobnum,k_sta:k_end)
complex(8) :: psi2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),  &
               1:iobnum,k_sta:k_end)
integer :: iter,p,q,iflag
integer :: ik
integer :: ix,iy,iz
integer :: is,pstart(2),pend(2)
complex(8) :: sum0,sum1,xkHxk,xkxk,Rk,gkgk,pkHpk
complex(8) :: uk
real(8) :: ev
complex(8) :: cx,cp,xkHpk,zs,xkTxk
!real(8) :: xk(ML),hxk(ML),gk(ML),pk(ML)
complex(8) , allocatable :: xk(:,:,:),hxk(:,:,:),gk(:,:,:),pk(:,:,:)
complex(8) , allocatable :: txk(:,:,:),htpsi(:,:,:),pko(:,:,:)
complex(8) , allocatable :: gk2(:,:,:)
real(8) :: elp2(2000)
complex(8):: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd)
integer :: p_myob,q_myob
integer :: icorr_p,icorr_q
integer :: iroot
integer :: is_sta,is_end

allocate (xk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (hxk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (gk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (gk2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (pk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

allocate (txk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (htpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (pko(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

call set_isstaend(is_sta,is_end)

iwk_size=2
call make_iwksta_iwkend

psi2=0.d0

!$OMP parallel do
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
  pstart(1)=1
  pend(1)=itotMST
else if(ilsda == 1)then
  pstart(1)=1
  pend(1)=MST(1)
  pstart(2)=MST(1)+1
  pend(2)=itotMST
end if

do ik=k_sta,k_end
do is=is_sta,is_end

orbital : do p=pstart(is),pend(is)
  call calc_myob(p,p_myob)
  call check_corrkob(p,ik,icorr_p)

  elp2(2)=get_wtime()

  if(nproc_ob==1)then
    do q=pstart(is),p-1
      sum0=0.d0
!$omp parallel do reduction(+ : sum0)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        sum0=sum0+conjg(psi_in(ix,iy,iz,q,ik))*psi_in(ix,iy,iz,p,ik)
      end do
      end do
      end do
      sum0=sum0*Hvol
      call comm_summation(sum0,sum1,nproc_group_korbital)
!$omp parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        psi_in(ix,iy,iz,p,ik)=psi_in(ix,iy,iz,p,ik)-sum1*psi_in(ix,iy,iz,q,ik)
      end do
      end do
      end do
    end do
  else
    do q=pstart(is),p-1
      call calc_myob(q,q_myob)
      call check_corrkob(q,ik,icorr_q)
      if(icorr_q==1)then
!$omp parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          cmatbox_m(ix,iy,iz)=psi_in(ix,iy,iz,q_myob,ik)
        end do
        end do
        end do
      end if
      call calc_iroot(q,iroot)
      call comm_bcast(cmatbox_m,nproc_group_kgrid,iroot)
      sum0=0.d0
      if(icorr_p==1)then
!$omp parallel do reduction(+ : sum0)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          sum0=sum0+conjg(cmatbox_m(ix,iy,iz))*psi_in(ix,iy,iz,p_myob,ik)
        end do
        end do
        end do
      end if
      sum0=sum0*Hvol
      call comm_summation(sum0,sum1,nproc_group_korbital)
      if(icorr_p==1)then
!$omp parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          psi_in(ix,iy,iz,p_myob,ik)=psi_in(ix,iy,iz,p_myob,ik)-sum1*cmatbox_m(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do 
  end if
  sum0=0.d0
  if(icorr_p==1)then
!$omp parallel do reduction(+ : sum0)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      sum0=sum0+abs(psi_in(ix,iy,iz,p_myob,ik))**2
    end do
    end do
    end do
  end if
  sum0=sum0*Hvol
  call comm_summation(sum0,sum1,nproc_group_korbital)
  if(icorr_p==1)then
!$omp parallel do 
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      xk(ix,iy,iz)=psi_in(ix,iy,iz,p_myob,ik)/sqrt(sum1)
      tpsi(ix,iy,iz)=xk(ix,iy,iz)
    end do
    end do
    end do
 
    call hpsi2(tpsi,hxk,p,ik,0,0)
  
!$omp parallel do 
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      psi2(ix,iy,iz,p_myob,ik)=hxk(ix,iy,iz)
    end do
    end do
    end do

!$omp parallel do 
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      txk(ix,iy,iz)=ttpsi(ix,iy,iz)
    end do
    end do
    end do
  end if

  call calc_iroot(p,iroot)
  call comm_bcast(xk,nproc_group_kgrid,iroot)
  call comm_bcast(hxk,nproc_group_kgrid,iroot)
  call comm_bcast(txk,nproc_group_kgrid,iroot)

  call inner_product4(xk,hxk,xkHxk)
  call inner_product4(xk,txk,xkTxk)
  
  Iteration : do iter=1,Ncg

!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      gk(ix,iy,iz) = hxk(ix,iy,iz) - xkHxk*xk(ix,iy,iz) 
    end do
    end do
    end do

    if(nproc_ob==1)then
      do q=pstart(is),p-1
        sum0=0.d0
!$omp parallel do reduction(+ : sum0)
        do iz=iwk3sta(3),iwk3end(3)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          sum0=sum0+conjg(psi_in(ix,iy,iz,q,ik))*gk(ix,iy,iz)
        end do
        end do
        end do
        sum0=sum0*Hvol
        call comm_summation(sum0,sum1,nproc_group_korbital)
        do iz=iwk3sta(3),iwk3end(3)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          gk(ix,iy,iz)=gk(ix,iy,iz)-sum1*psi_in(ix,iy,iz,q,ik)
        end do
        end do
        end do
      end do
    else
      do q=pstart(is),p-1
        call calc_myob(q,q_myob)
        call check_corrkob(q,ik,icorr_q)
        if(icorr_q==1)then
!$omp parallel do
          do iz=mg_sta(3),mg_end(3)
          do iy=mg_sta(2),mg_end(2)
          do ix=mg_sta(1),mg_end(1)
            cmatbox_m(ix,iy,iz)=psi_in(ix,iy,iz,q_myob,ik)
          end do
          end do
          end do
        end if
        call calc_iroot(q,iroot)
        call comm_bcast(cmatbox_m,nproc_group_kgrid,iroot)
        sum0=0.d0
!$omp parallel do reduction(+ : sum0)
        do iz=iwk3sta(3),iwk3end(3)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          sum0=sum0+conjg(cmatbox_m(ix,iy,iz))*gk(ix,iy,iz)
        end do
        end do
        end do
        sum0=sum0*Hvol
        call comm_summation(sum0,sum1,nproc_group_korbital)
        do iz=iwk3sta(3),iwk3end(3)
        do iy=iwk3sta(2),iwk3end(2)
        do ix=iwk3sta(1),iwk3end(1)
          gk(ix,iy,iz)=gk(ix,iy,iz)-sum1*cmatbox_m(ix,iy,iz)
        end do
        end do
        end do
      end do
    end if
    call inner_product4(gk,gk,sum1)
    
    if(iter==1)then
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        pk(ix,iy,iz)=gk(ix,iy,iz)
      end do
      end do
      end do
    else
      uk=sum1/gkgk
!$OMP parallel do
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        pk(ix,iy,iz)=gk(ix,iy,iz)+uk*pk(ix,iy,iz)
      end do
      end do
      end do
    end if
    gkgk=sum1
    call inner_product4(xk,pk,zs)
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      pko(ix,iy,iz)=pk(ix,iy,iz)-zs*xk(ix,iy,iz)
    end do
    end do
    end do
    call inner_product4(pko,pko,sum1)
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      pko(ix,iy,iz)=pko(ix,iy,iz)/sqrt(sum1)
      tpsi(ix,iy,iz)=pko(ix,iy,iz)
    end do
    end do
    end do
    call hpsi2(tpsi,htpsi,p,ik,0,0)
    call inner_product4(xk,htpsi,xkHpk)
    call inner_product4(pko,htpsi,pkHpk)
    

    ev=0.5d0*((xkHxk+pkHpk)-sqrt((xkHxk-pkHpk)**2+4.d0*abs(xkHpk)**2))
    cx=xkHpk/(ev-xkHxk)
    cp=1.d0/sqrt(1.d0+abs(cx)**2)
    cx=cx*cp

    
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      xk(ix,iy,iz)=cx*xk(ix,iy,iz)+cp*pko(ix,iy,iz)
      hxk(ix,iy,iz)=cx*hxk(ix,iy,iz)+cp*htpsi(ix,iy,iz)
      txk(ix,iy,iz)=cx*txk(ix,iy,iz)+cp*ttpsi(ix,iy,iz)
    end do
    end do
    end do

    call inner_product4(xk,hxk,xkHxk)
    call inner_product4(xk,txk,xkTxk)
    call inner_product4(xk,xk,xkxk)
    Rk=xkHxk/xkxk

  end do Iteration

  call inner_product4(xk(mg_sta(1),mg_sta(2),mg_sta(3)),xk(mg_sta(1),mg_sta(2),mg_sta(3)),sum0)
  if(icorr_p==1)then
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      psi_in(ix,iy,iz,p_myob,ik)=xk(ix,iy,iz)/sqrt(sum0)
    end do
    end do
    end do
  end if
end do orbital

end do
end do

if(iflag.eq.1) then
  iflag=0
end if

deallocate (xk,hxk,gk,pk,gk2)

return

END SUBROUTINE DTcg_periodic

