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

SUBROUTINE DTcg(psi_in,iflag)
use salmon_parallel, only: nproc_group_grid
use salmon_communication, only: comm_bcast
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
real(8) :: sum0,xkHxk,xkxk,Rk,gkgk,xkpk,pkpk,pkHxk,pkHpk
real(8) :: uk,alpha,Ak,Bk,Ck
real(8) , allocatable :: xk(:,:,:),hxk(:,:,:),gk(:,:,:),pk(:,:,:)
real(8) , allocatable :: gk2(:,:,:)
real(8) :: elp2(2000)
real(8):: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd)
integer :: iob_myob,job_myob
integer :: icorr,jcorr               
integer :: iroot
integer :: is_sta,is_end

allocate (xk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (hxk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (gk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (gk2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate (pk(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

iwk_size=2
call make_iwksta_iwkend

call set_isstaend(is_sta,is_end)

!$OMP parallel do private(iz,iy,ix) collapse(3)
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

do is=is_sta,is_end


orbital : do iob=iobsta(is),iobend(is)
  call calc_myob(iob,iob_myob)
  call check_corrkob(iob,icorr)
  elp2(2)=get_wtime()

  if(icorr==1)then

!$OMP parallel do private(iz,iy,ix) collapse(3)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      xk(ix,iy,iz)=psi_in(ix,iy,iz,iob_myob,1)
    end do
    end do
    end do

!$OMP parallel do private(iz,iy,ix) collapse(3)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      tpsi(ix,iy,iz)=xk(ix,iy,iz)
    end do
    end do
    end do

    call hpsi2(tpsi,hxk,iob,0,0)

    call inner_product(xk,hxk,xkHxk)

    xkHxk=xkHxk*Hvol ; xkxk=1.d0 ; Rk=xkHxk/xkxk

  end if

  Iteration : do iter=1,Ncg

    if(icorr==1)then
!$OMP parallel do private(iz,iy,ix) collapse(3)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        gk(ix,iy,iz) = 2*( hxk(ix,iy,iz) - Rk*xk(ix,iy,iz) )
      end do
      end do
      end do
    end if
    call calc_iroot(iob,iroot)
    call comm_bcast(gk,nproc_group_grid,iroot)

    do job=iobsta(is),iob-1
      sum0=0.d0
      call calc_myob(job,job_myob)
      call check_corrkob(job,jcorr)
      if(jcorr==1)then
        call inner_product(psi_in(:,:,:,job_myob,1),gk(:,:,:),sum0)
        sum0=sum0*Hvol
!$OMP parallel do private(iz,iy,ix) collapse(3)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          gk(ix,iy,iz)=gk(ix,iy,iz)-sum0*psi_in(ix,iy,iz,job_myob,1)
        end do
        end do
        end do
      end if

      call calc_iroot(job,iroot)
      call comm_bcast(gk,nproc_group_grid,iroot)
    end do

    if(icorr==1)then

      call inner_product(gk,gk,sum0)
      sum0=sum0*Hvol

      if ( iter==1 ) then
!$OMP parallel do private(iz,iy,ix) collapse(3)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          pk(ix,iy,iz) = -gk(ix,iy,iz)
        end do
        end do
        end do
      else
        uk=sum0/gkgk 
!$OMP parallel do private(iz,iy,ix) collapse(3)
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          pk(ix,iy,iz) = -gk(ix,iy,iz) + uk*pk(ix,iy,iz)
        end do
        end do
        end do
      end if
      gkgk=sum0

      xkpk=0.d0 ; pkpk=0.d0 ; pkHxk=0.d0

      call inner_product(xk,pk,xkpk)
      xkpk = xkpk*Hvol

      call inner_product(pk,pk,pkpk)
      pkpk = pkpk*Hvol

      call inner_product(pk,hxk,pkHxk)
      pkHxk = pkHxk*Hvol

!$OMP parallel do private(iz,iy,ix) collapse(3)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        tpsi(ix,iy,iz) = pk(ix,iy,iz)
      end do
      end do
      end do
      call hpsi2(tpsi,gk,iob,0,0)

      call inner_product(pk,gk,pkHpk)
      pkHpk = pkHpk*Hvol

      Ak=pkHpk*xkpk-pkHxk*pkpk
      Bk=pkHpk*xkxk-xkHxk*pkpk
      Ck=pkHxk*xkxk-xkHxk*xkpk
      alpha=(-Bk+sqrt(Bk*Bk-4*Ak*Ck))/(2*Ak)

      xk = xk + alpha*pk
      hxk=hxk + alpha*gk

      call inner_product(xk,hxk,xkHxk)
      xkHxk = xkHxk*Hvol

      call inner_product(xk,xk,xkxk)
      xkxk = xkxk*Hvol
    
      Rk=xkHxk/xkxk

    end if

  end do Iteration

  if(icorr==1)then
    call inner_product(xk,xk,sum0)
!$OMP parallel do private(iz,iy,ix) collapse(3)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      psi_in(ix,iy,iz,iob_myob,1)=xk(ix,iy,iz)/sqrt(sum0*Hvol)
    end do
    end do
    end do
  end if

end do orbital

end do

if(iflag.eq.1) then
   iflag=0
end if

deallocate (xk,hxk,gk,pk,gk2)
return

END SUBROUTINE DTcg
