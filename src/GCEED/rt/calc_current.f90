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
subroutine calc_current(tpsi)
use salmon_parallel, only: nproc_group_global, nproc_group_korbital
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use scf_data
use sendrecv_groupob_tmp_sub
use allocate_psl_sub
implicit none
complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd, &
                   1:iobnum,k_sta:k_end)
integer :: ix,iy,iz,iob,iik
complex(8),parameter :: zi=(0.d0,1.d0)
complex(8) :: ekr(maxMps,MI,k_sta:k_end)
real(8) :: r(3)
integer :: iatom,ik,jj,lm
complex(8) :: uVpsi0,uVpsix,uVpsiy,uVpsiz
real(8) :: jxt,jyt,jzt
real(8) :: curr1(3),curr2(3)
integer :: p_allob

iwk_size=2
call make_iwksta_iwkend

curr1(1:3)=0.d0

elp3(1301)=get_wtime()

call sendrecv_groupob_tmp(tpsi)

elp3(1302)=get_wtime()
elp3(1351)=elp3(1351)+elp3(1302)-elp3(1301)

do iik=k_sta,k_end
do iob=1,iobnum
  call calc_allob(iob,p_allob)
  jxt=0.d0
!$OMP parallel do private(iz,iy,ix) reduction(+ : jxt)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    jxt=jxt+real(tpsi(ix,iy,iz,iob,iik))**2+aimag(tpsi(ix,iy,iz,iob,iik))**2
  end do
  end do
  end do
  jyt=jxt ; jzt=jxt
  curr1(1)=curr1(1)+rocc(p_allob,iik)*wtk(iik)*k_rd(1,iik)*jxt
  curr1(2)=curr1(2)+rocc(p_allob,iik)*wtk(iik)*k_rd(2,iik)*jyt
  curr1(3)=curr1(3)+rocc(p_allob,iik)*wtk(iik)*k_rd(3,iik)*jzt

  jxt=0.d0; jyt=0.d0; jzt=0.d0
!$OMP parallel do private(iz,iy,ix) reduction(+ : jxt,jyt,jzt)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    jxt=jxt+aimag(conjg(tpsi(ix,iy,iz,iob,iik))*   &
                            ( bN1*( tpsi(ix+1,iy,iz,iob,iik) - tpsi(ix-1,iy,iz,iob,iik) )    &
                             +bN2*( tpsi(ix+2,iy,iz,iob,iik) - tpsi(ix-2,iy,iz,iob,iik) )    &
                             +bN3*( tpsi(ix+3,iy,iz,iob,iik) - tpsi(ix-3,iy,iz,iob,iik) )    &
                             +bN4*( tpsi(ix+4,iy,iz,iob,iik) - tpsi(ix-4,iy,iz,iob,iik) ))/Hgs(1))
    jyt=jyt+aimag(conjg(tpsi(ix,iy,iz,iob,iik))*   &
                            ( bN1*( tpsi(ix,iy+1,iz,iob,iik) - tpsi(ix,iy-1,iz,iob,iik) )    &
                             +bN2*( tpsi(ix,iy+2,iz,iob,iik) - tpsi(ix,iy-2,iz,iob,iik) )    &
                             +bN3*( tpsi(ix,iy+3,iz,iob,iik) - tpsi(ix,iy-3,iz,iob,iik) )    &
                             +bN4*( tpsi(ix,iy+4,iz,iob,iik) - tpsi(ix,iy-4,iz,iob,iik) ))/Hgs(2))
    jzt=jzt+aimag(conjg(tpsi(ix,iy,iz,iob,iik))*   &
                            ( bN1*( tpsi(ix,iy,iz+1,iob,iik) - tpsi(ix,iy,iz-1,iob,iik) )    &
                             +bN2*( tpsi(ix,iy,iz+2,iob,iik) - tpsi(ix,iy,iz-2,iob,iik) )    &
                             +bN3*( tpsi(ix,iy,iz+3,iob,iik) - tpsi(ix,iy,iz-3,iob,iik) )    &
                             +bN4*( tpsi(ix,iy,iz+4,iob,iik) - tpsi(ix,iy,iz-4,iob,iik) ))/Hgs(3))
  end do
  end do
  end do
  curr1(1)=curr1(1)+rocc(p_allob,iik)*wtk(iik)*jxt
  curr1(2)=curr1(2)+rocc(p_allob,iik)*wtk(iik)*jyt
  curr1(3)=curr1(3)+rocc(p_allob,iik)*wtk(iik)*jzt
end do
end do
curr1(1:3)=curr1(1:3)*Hvol/(dble(lg_num(1)*lg_num(2)*lg_num(3))*Hvol)

elp3(1303)=get_wtime()
elp3(1352)=elp3(1352)+elp3(1303)-elp3(1302)

jxt=0.d0;jyt=0.d0;jzt=0.d0
do iik=k_sta,k_end
!$OMP parallel do private(iatom,jj,ik,r)
  do iatom=1,MI
    ik=Kion(iatom)
    do jj=1,Mps(iatom)
      r(1)=(dble(Jxyz(1,jj,iatom)-1)-dble(Jxxyyzz(1,jj,iatom)*lg_num(1)))*Hgs(1)
      r(2)=(dble(Jxyz(2,jj,iatom)-1)-dble(Jxxyyzz(2,jj,iatom)*lg_num(2)))*Hgs(2)
      r(3)=(dble(Jxyz(3,jj,iatom)-1)-dble(Jxxyyzz(3,jj,iatom)*lg_num(3)))*Hgs(3)
      ekr(jj,iatom,iik)=exp(zi*(k_rd(1,iik)*r(1)   &
                          +k_rd(2,iik)*r(2)   &
                          +k_rd(3,iik)*r(3)))
    end do
  end do
  do iob=1,iobnum
    call calc_allob(iob,p_allob)
!$OMP parallel do private(iatom,lm,ik,uVpsi0,uVpsix,uVpsiy,uVpsiz,r)
    do iatom=1,MI
      ik=Kion(iatom)
      do lm=1,(Mlps(ik)+1)**2
        if ( abs(uVu(lm,iatom))<1.d-5 ) then
          uVpsibox1_j(1:4,lm,iatom,iob,iik)=0.d0
        else
          uVpsi0=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
          do jj=1,Mps(iatom)
            r(1)=(dble(Jxyz(1,jj,iatom)-1)-dble(Jxxyyzz(1,jj,iatom)*lg_num(1)))*Hgs(1)
            r(2)=(dble(Jxyz(2,jj,iatom)-1)-dble(Jxxyyzz(2,jj,iatom)*lg_num(2)))*Hgs(2)
            r(3)=(dble(Jxyz(3,jj,iatom)-1)-dble(Jxxyyzz(3,jj,iatom)*lg_num(3)))*Hgs(3)
            uVpsi0=uVpsi0+uV(jj,lm,iatom)*ekr(jj,iatom,iik)      &
                     *tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
            uVpsix=uVpsix+uV(jj,lm,iatom)*ekr(jj,iatom,iik)*r(1) &
                     *tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
            uVpsiy=uVpsiy+uV(jj,lm,iatom)*ekr(jj,iatom,iik)*r(2) &
                     *tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
            uVpsiz=uVpsiz+uV(jj,lm,iatom)*ekr(jj,iatom,iik)*r(3) &
                     *tpsi(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)
          end do  
          uVpsibox1_j(1,lm,iatom,iob,iik)=uVpsi0*Hvol/uVu(lm,iatom)
          uVpsibox1_j(2,lm,iatom,iob,iik)=uVpsix*Hvol
          uVpsibox1_j(3,lm,iatom,iob,iik)=uVpsiy*Hvol
          uVpsibox1_j(4,lm,iatom,iob,iik)=uVpsiz*Hvol
        end if
      end do
    end do
  end do
end do
elp3(1304)=get_wtime()
elp3(1353)=elp3(1353)+elp3(1304)-elp3(1303)

call comm_summation(uVpsibox1_j,uVpsibox2_j,4*maxlm*MI*iobnum*k_num,nproc_group_korbital)
elp3(1305)=get_wtime()
elp3(1354)=elp3(1354)+elp3(1305)-elp3(1304)

do iik=k_sta,k_end
  do iob=1,iobnum
    call calc_allob(iob,p_allob)
    do iatom=1,MI
      ik=Kion(iatom)
      do lm=1,(Mlps(ik)+1)**2
        jxt=jxt+rocc(p_allob,iik)*wtk(iik)/(lg_num(1)*lg_num(2)*lg_num(3)*Hvol)   &
                          *2.d0*aimag(conjg(uVpsibox2_j(2,lm,iatom,iob,iik))  &
                          *uVpsibox2_j(1,lm,iatom,iob,iik))/dble(nproc_Mxin_mul) 
        jyt=jyt+rocc(p_allob,iik)*wtk(iik)/(lg_num(1)*lg_num(2)*lg_num(3)*Hvol)   &
                          *2.d0*aimag(conjg(uVpsibox2_j(3,lm,iatom,iob,iik))  &
                          *uVpsibox2_j(1,lm,iatom,iob,iik))/dble(nproc_Mxin_mul)
        jzt=jzt+rocc(p_allob,iik)*wtk(iik)/(lg_num(1)*lg_num(2)*lg_num(3)*Hvol)   &
                          *2.d0*aimag(conjg(uVpsibox2_j(4,lm,iatom,iob,iik))  &
                          *uVpsibox2_j(1,lm,iatom,iob,iik))/dble(nproc_Mxin_mul)
      end do
    end do
  end do
end do

curr1(1)=curr1(1)+jxt
curr1(2)=curr1(2)+jyt
curr1(3)=curr1(3)+jzt

elp3(1306)=get_wtime()
elp3(1355)=elp3(1355)+elp3(1306)-elp3(1305)

call comm_summation(curr1,curr2,3,nproc_group_global)

elp3(1307)=get_wtime()
elp3(1356)=elp3(1356)+elp3(1307)-elp3(1306)

curr(1:3,itt)=curr2(1:3)

if(iflag_indA==1)then
  A_ind(:,itt+1)=2.d0*A_ind(:,itt)-A_ind(:,itt-1)-4.d0*Pi*curr(:,itt)*dt**2
else if(iflag_indA==0)then
  A_ind(:,itt+1)=0.d0
end if

A_tot(:,itt+1)=A_ext(:,itt+1)+A_ind(:,itt+1)

E_ext(:,itt)=-(A_ext(:,itt+1)-A_ext(:,itt-1))/(2.d0*dt)
E_ind(:,itt)=-(A_ind(:,itt+1)-A_ind(:,itt-1))/(2.d0*dt)
E_tot(:,itt)=-(A_tot(:,itt+1)-A_tot(:,itt-1))/(2.d0*dt)

end subroutine calc_current
