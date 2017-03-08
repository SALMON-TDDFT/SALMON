!
!  Copyright 2016 ARTED developers
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
!This file is "CG.f90"
!This file contain a subroutine.
!Subroutine CG
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine CG_omp(iter_cg_max)
  use Global_Variables
  use communication
  use timelog
  implicit none
  real(8),parameter :: delta_cg=1.d-15
  integer iter,ik,ib,ibt
  integer :: iter_cg_max
  real(8) :: xkHxk,gkgk,pkHpk,xkTxk
  real(8) :: uk,s,ev
!  complex(8) :: xk(NL),hxk(NL),gk(NL),pk(NL),pko(NL),txk(NL)
  complex(8) :: cx,cp,xkHpk
  real(8) :: esp_var_l(1:NB,1:NK)
  complex(8) :: zs
! sato
  integer :: ia,j,i,ix,iy,iz
  real(8) :: kr
! omp
  integer :: thr_id,omp_get_thread_num
  thr_id=0

  call timelog_begin(LOG_CG)
  esp_var_l(:,:)=0.d0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ia,j,i,ix,iy,iz,kr,ib,ibt,s,xkHxk,xkTxk,iter,uk,gkgk,xkHpk,pkHpk,ev,cx,cp,zs)
  do ik=NK_s,NK_e
! sato -----------------------------------------------------------------------------------------
!    k2=sum(kAc(ik,:)**2)
!Constructing nonlocal part
  do ia=1,NI
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  enddo
  enddo
! sato -----------------------------------------------------------------------------------------
  do ib=1,NB
    do ibt=1,ib-1
      s=sum(conjg(zu_GS(:,ibt,ik))*zu_GS(:,ib,ik))*Hxyz
      zu_GS(1:NL,ib,ik)=zu_GS(1:NL,ib,ik)-zu_GS(1:NL,ibt,ik)*s
    end do
    s=1.0d0/sqrt(sum(abs(zu_GS(:,ib,ik))**2)*Hxyz)
    xk_omp(1:NL,thr_id)=zu_GS(1:NL,ib,ik)*s
    call hpsi_omp_KB_GS(ik,xk_omp(:,thr_id),txk_omp(:,thr_id),hxk_omp(:,thr_id))
    xkHxk=sum(conjg(xk_omp(:,thr_id))*hxk_omp(:,thr_id))*Hxyz
    xkTxk=sum(conjg(xk_omp(:,thr_id))*txk_omp(:,thr_id))*Hxyz

    do iter=1,iter_cg_max
      gk_omp(1:NL,thr_id)=(hxk_omp(1:NL,thr_id)-xkHxk*xk_omp(1:NL,thr_id))
      do ibt=1,ib-1
        zs=sum(conjg(zu_GS(:,ibt,ik))*gk_omp(:,thr_id))*Hxyz
        gk_omp(1:NL,thr_id)=gk_omp(1:NL,thr_id)-zu_GS(1:NL,ibt,ik)*zs
      end do
      s=sum(abs(gk_omp(:,thr_id))**2)*Hxyz

      select case (iter)
      case(1)
        pk_omp(1:NL,thr_id)=gk_omp(1:NL,thr_id)
      case default
        uk=s/gkgk
        pk_omp(1:NL,thr_id)=gk_omp(1:NL,thr_id)+uk*pk_omp(1:NL,thr_id)
      end select
      gkgk=s

      zs=sum(conjg(xk_omp(:,thr_id))*pk_omp(:,thr_id))*Hxyz
      pko_omp(1:NL,thr_id)=pk_omp(1:NL,thr_id)-xk_omp(1:NL,thr_id)*zs
      s=1.0d0/sqrt(sum(abs(pko_omp(:,thr_id))**2)*Hxyz)
      pko_omp(1:NL,thr_id)=pko_omp(1:NL,thr_id)*s
      call hpsi_omp_KB_GS(ik,pko_omp(:,thr_id),ttpsi_omp(:,thr_id),htpsi_omp(:,thr_id))
      xkHpk=sum(conjg( xk_omp(:,thr_id))*htpsi_omp(:,thr_id))*Hxyz
      pkHpk=sum(conjg(pko_omp(:,thr_id))*htpsi_omp(:,thr_id))*Hxyz
      ev=0.5d0*((xkHxk+pkHpk)-sqrt((xkHxk-pkHpk)**2+4*abs(xkHpk)**2))
      cx=xkHpk/(ev-xkHxk)
      cp=1.d0/sqrt(1.d0+abs(cx)**2)
      cx=cx*cp
      if(abs(ev-xkHxk)<delta_cg) exit
       xk_omp(1:NL,thr_id)=cx* xk_omp(1:NL,thr_id)+cp*  pko_omp(1:NL,thr_id)
      hxk_omp(1:NL,thr_id)=cx*hxk_omp(1:NL,thr_id)+cp*htpsi_omp(1:NL,thr_id)
      txk_omp(1:NL,thr_id)=cx*txk_omp(1:NL,thr_id)+cp*ttpsi_omp(1:NL,thr_id)
      xkHxk=sum(conjg(xk_omp(:,thr_id))*hxk_omp(:,thr_id))*Hxyz
      xkTxk=sum(conjg(xk_omp(:,thr_id))*txk_omp(:,thr_id))*Hxyz
    enddo

    s=1.0d0/sqrt(sum(abs(xk_omp(:,thr_id))**2)*Hxyz)
    zu_GS(1:NL,ib,ik)=xk_omp(1:NL,thr_id)*s
    call hpsi_omp_KB_GS(ik,zu_GS(:,ib,ik),ttpsi_omp(:,thr_id),htpsi_omp(:,thr_id))
    xkHxk=sum(conjg(zu_GS(1:NL,ib,ik))*htpsi_omp(:,thr_id))*Hxyz
    esp_var_l(ib,ik)=sqrt(sum(abs(htpsi_omp(:,thr_id)-xkHxk*zu_GS(1:NL,ib,ik))**2)*Hxyz)*occ(ib,ik)
  enddo
  enddo

!$omp end parallel

  call timelog_begin(LOG_ALLREDUCE)
  call comm_summation(esp_var_l,esp_var,NB*NK,proc_group(2))
  call timelog_end(LOG_ALLREDUCE)

  call timelog_end(LOG_CG)

  return
End Subroutine CG_OMP
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
