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
subroutine Total_Energy_periodic_scf(tzpsi_in)
use salmon_parallel, only: nproc_group_global, nproc_id_global, nproc_size_global, nproc_group_rho, nproc_group_korbital
use salmon_communication, only: comm_is_root, comm_summation
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
use read_pslfile_sub
use allocate_psl_sub
use allocate_mat_sub
use sendrecv_tmp_sub
use sendrecv_groupob_ngp_sub
implicit none

complex(8) :: tzpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd,  &
                1:iobnum,k_sta:k_end)

complex(8),parameter :: zi=(0.d0,1.d0)
integer :: ii,jj,iatom,iob,ia,ib,lm,ind
integer :: ia_sta,ia_end
real(8) :: x,y,z
integer :: ik,iik
integer :: ix,iy,iz
integer :: iix,iiy,iiz
integer :: aewald_sta,aewald_end,aewald_num
integer :: totnum_aewald
real(8) :: rab(3),rab2
real(8) :: rbox1,rbox2,rbox3,rbox4
real(8) :: Ekin_tmp
real(8) :: Ebox1(9),Ebox2(9) ! 1:Ekin, 2:Eion, 3:Eion_l, 4:Eh_l, 5:Eloc_l1, 6: Eloc_l2, 7:Enl_l
real(8) :: bLx,bLy,bLz
integer :: NG_s,NG_e
integer :: NG_l_s_para,NG_l_e_para
integer :: numtmp
real(8) :: G2,Gd
real(8) :: sysvol
integer(8) :: n
complex(8) :: uVpsibox3(1:maxlm,1:MI,1:iobnum,k_sta:k_end)
complex(8) :: uVpsibox4(1:maxlm,1:MI,1:iobnum,k_sta:k_end)
complex(8) :: ekr(maxMps,MI,k_sta:k_end)
complex(8) :: sumbox

complex(8) :: fdN0
complex(8) :: fdN1(0:12,3)
complex(8) :: fdN2(0:12,3)
real(8) :: f0
integer :: p_allob

complex(8) :: htpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
complex(8) :: tpsi_tmp2(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd)

real(8) :: sum_temp1,sum_temp2
real(8) :: sum_temp3(4)
real(8) :: sum_temp4(4)

integer :: iy_sta,iy_end,iz_sta,iz_end
real(8) :: rbox13,rbox14,rbox15,rbox16

!if(t>=11) call fapp_start("region8",8,0)

iwk_size=2
call make_iwksta_iwkend

if(iperiodic==3)then

f0=(1.d0/Hgs(1)**2   &
   +1.d0/Hgs(2)**2   &
   +1.d0/Hgs(3)**2)

elp3(1401)=get_wtime()

tpsi_tmp2=0.d0
do iik=k_sta,k_end
do iob=1,iobnum
!$OMP parallel do private(ix,iy,iz)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    tpsi_tmp2(ix,iy,iz)= tzpsi_in(ix,iy,iz,iob,iik)
  end do
  end do
  end do
  call sendrecv_tmp(tpsi_tmp2)
!$OMP parallel do private(ix,iy,iz)
  do iz=mg_sta(3)-Nd,mg_end(3)+Nd
  do iy=mg_sta(2)-Nd,mg_end(2)+Nd
  do ix=mg_sta(1)-Nd,mg_end(1)+Nd
    tzpsi_in(ix,iy,iz,iob,iik) = tpsi_tmp2(ix,iy,iz)
  end do
  end do
  end do
end do
end do

elp3(1402)=get_wtime()
elp3(1451)=elp3(1451)+elp3(1402)-elp3(1401)

Ebox1(1:9)=0.d0

do iik=k_sta,k_end
  fdN0=-0.5d0*cNmat(0,Nd)*f0
  do jj=1,3
    do ind=1,4
       fdN1(ind,jj)=-0.5d0*cNmat(ind,Nd)/Hgs(jj)**2-zi*k_rd(jj,iik)*bNmat(ind,Nd)/Hgs(jj)
       fdN2(ind,jj)=-0.5d0*cNmat(ind,Nd)/Hgs(jj)**2+zi*k_rd(jj,iik)*bNmat(ind,Nd)/Hgs(jj)
    end do
  end do
  do iob=1,iobnum
    call calc_allob(iob,p_allob)
!$OMP parallel do private(ix,iy,iz)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      htpsi(ix,iy,iz) =  &
         fdN0*tzpsi_in(ix,iy,iz,iob,iik)      &
        +fdN1(1,1)* tzpsi_in(ix+1,iy,iz,iob,iik) + fdN2(1,1) * tzpsi_in(ix-1,iy,iz,iob,iik)      &
        +fdN1(2,1)* tzpsi_in(ix+2,iy,iz,iob,iik) + fdN2(2,1) * tzpsi_in(ix-2,iy,iz,iob,iik)      &
        +fdN1(3,1)* tzpsi_in(ix+3,iy,iz,iob,iik) + fdN2(3,1) * tzpsi_in(ix-3,iy,iz,iob,iik)      &
        +fdN1(4,1)* tzpsi_in(ix+4,iy,iz,iob,iik) + fdN2(4,1) * tzpsi_in(ix-4,iy,iz,iob,iik)      &
        +fdN1(1,2)* tzpsi_in(ix,iy+1,iz,iob,iik) + fdN2(1,2) * tzpsi_in(ix,iy-1,iz,iob,iik)      &
        +fdN1(2,2)* tzpsi_in(ix,iy+2,iz,iob,iik) + fdN2(2,2) * tzpsi_in(ix,iy-2,iz,iob,iik)      &
        +fdN1(3,2)* tzpsi_in(ix,iy+3,iz,iob,iik) + fdN2(3,2) * tzpsi_in(ix,iy-3,iz,iob,iik)      &
        +fdN1(4,2)* tzpsi_in(ix,iy+4,iz,iob,iik) + fdN2(4,2) * tzpsi_in(ix,iy-4,iz,iob,iik)      &
        +fdN1(1,3)* tzpsi_in(ix,iy,iz+1,iob,iik) + fdN2(1,3) * tzpsi_in(ix,iy,iz-1,iob,iik)      &
        +fdN1(2,3)* tzpsi_in(ix,iy,iz+2,iob,iik) + fdN2(2,3) * tzpsi_in(ix,iy,iz-2,iob,iik)      &
        +fdN1(3,3)* tzpsi_in(ix,iy,iz+3,iob,iik) + fdN2(3,3) * tzpsi_in(ix,iy,iz-3,iob,iik)      &
        +fdN1(4,3)* tzpsi_in(ix,iy,iz+4,iob,iik) + fdN2(4,3) * tzpsi_in(ix,iy,iz-4,iob,iik) 
    end do
    end do
    end do
    Ekin_tmp=0.d0
!$OMP parallel do reduction(+:Ekin_tmp) private(ix,iy,iz)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Ekin_tmp=Ekin_tmp + rocc(iob,iik)*wtk(iik)*conjg(tzpsi_in(ix,iy,iz,iob,iik))*htpsi(ix,iy,iz)*Hvol 
    end do
    end do
    end do
    Ebox1(1)=Ebox1(1)+Ekin_tmp + rocc(p_allob,iik)*wtk(iik)*(ksquare(iik))/2.d0/dble(nproc_Mxin_mul)
  end do
end do

sum_temp1=Ebox1(1)
call comm_summation(sum_temp1,sum_temp2,nproc_group_global)
Ebox2(1)=sum_temp2

elp3(1403)=get_wtime()
elp3(1452)=elp3(1452)+elp3(1403)-elp3(1402)

totnum_aewald=MI*(NEwald*2+1)**3

aewald_sta=nproc_id_global*totnum_aewald/nproc_size_global+1
aewald_end=(nproc_id_global+1)*totnum_aewald/nproc_size_global
aewald_num=aewald_end-aewald_sta+1

!    Eion_tmp1=Eion_tmp1-Pi*sum(Zps(Kion(:)))**2/(2*aEwald*aLxyz) - sqrt(aEwald/Pi)*sum(Zps(Kion(:))**2)
if(aewald_num>=1)then
! 2:Eion
  do ia=1,MI
  do ii=aewald_sta,aewald_end
    iix=mod((ii-1)/((NEwald*2+1)**2),NEwald*2+1)-NEwald
    iiy=mod((ii-1)/(NEwald*2+1),NEwald*2+1)-NEwald
    iiz=mod(ii-1,NEwald*2+1)-NEwald
    do ib=1,MI
      if (iix**2+iiy**2+iiz**2 == 0 .and. ia == ib) cycle
      rab(1)=Rion(1,ia)-dble(iix)*Hgs(1)*lg_num(1)-Rion(1,ib)
      rab(2)=Rion(2,ia)-dble(iiy)*Hgs(2)*lg_num(2)-Rion(2,ib)
      rab(3)=Rion(3,ia)-dble(iiz)*Hgs(3)*lg_num(3)-Rion(3,ib)
      rab2=sum(rab(:)**2)
      Ebox1(2)=Ebox1(2) + 0.5d0*Zps(Kion(ia))*Zps(Kion(ib))*erfc(sqrt(aEwald*rab2))/sqrt(rab2)
    enddo
  enddo
  enddo
end if

sum_temp1=Ebox1(2)
call comm_summation(sum_temp1,sum_temp2,nproc_group_global)
Ebox2(2)=sum_temp2

ia_sta=nproc_id_global*MI/nproc_size_global+1
ia_end=(nproc_id_global+1)*MI/nproc_size_global
rbox1=0.d0
rbox2=0.d0
do ia=ia_sta,ia_end
  ik=Kion(ia)
  rbox1=rbox1+Zps(ik)
  rbox2=rbox2+Zps(ik)**2
end do
call comm_summation(rbox1,rbox3,nproc_group_global)
call comm_summation(rbox2,rbox4,nproc_group_global)
Ebox2(2)=Ebox2(2)-Pi*rbox3**2/(2*aEwald*Hvol*lg_num(1)*lg_num(2)*lg_num(3))  &
                 - sqrt(aEwald/Pi)*rbox4


elp3(1404)=get_wtime()
elp3(1453)=elp3(1453)+elp3(1404)-elp3(1403)

!calculate reciprocal lattice vector
bLx=2.d0*Pi/(Hgs(1)*dble(lg_num(1)))
bLy=2.d0*Pi/(Hgs(2)*dble(lg_num(2)))
bLz=2.d0*Pi/(Hgs(3)*dble(lg_num(3)))

NG_s=1
NG_e=lg_num(1)*lg_num(2)*lg_num(3)

numtmp=NG_e/nproc_size_global
NG_l_s_para = nproc_id_global*numtmp+1
NG_l_e_para = (nproc_id_global+1)*numtmp
if(nproc_id_global==nproc_size_global-1) NG_l_e_para=NG_e

sysvol=Hvol*lg_num(1)*lg_num(2)*lg_num(3)

!3:Eion_l, 4:Eh_l, 5:Eloc_l1, 6: Eloc_l2
select case(iflag_hartree)
case(2)
  do n=NG_l_s_para,NG_l_e_para
    if(n == nGzero ) cycle
    G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
    Ebox1(3)=Ebox1(3)+sysvol*(4*Pi/G2)*(abs(rhoion_G(n))**2*exp(-G2/(4*aEwald))*0.5d0)
    Ebox1(4)=Ebox1(4)+sysvol*(4*Pi/G2)*(abs(rhoe_G(n))**2*0.5d0)
    Ebox1(5)=Ebox1(5)+sysvol*(4*Pi/G2)*(-rhoe_G(n)*conjg(rhoion_G(n)))
    do ia=1,MI
      Gd=Gx(n)*Rion(1,ia)+Gy(n)*Rion(2,ia)+Gz(n)*Rion(3,ia)
      Ebox1(6)=Ebox1(6)+conjg(rhoe_G(n))*dVloc_G(n,Kion(ia))*exp(-zI*Gd)
    end do
  enddo
case(4)
  iz_sta=1
  iz_end=lg_num(3)/NPUZ
  iy_sta=1
  iy_end=lg_num(2)/NPUY

  do iz=iz_sta,iz_end
    do iy=iy_sta,iy_end
      rbox13=0.d0
      rbox14=0.d0
      rbox15=0.d0
      rbox16=0.d0
!$OMP parallel do reduction (+ : rbox13,rbox14,rbox15,rbox16) private(n,G2,iix,Gd)
      do ix=1,lg_num(1)
        n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
        if(n == nGzero ) cycle
        G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
        rbox13=rbox13+sysvol*(4*Pi/G2)*(abs(rhoion_G(n))**2*exp(-G2/(4*aEwald))*0.5d0)
        rbox14=rbox14+sysvol*(4*Pi/G2)*(abs(rhoe_G(n))**2*0.5d0)
        rbox15=rbox15+sysvol*(4*Pi/G2)*(-rhoe_G(n)*conjg(rhoion_G(n)))
        do ia=1,MI
          Gd=Gx(n)*Rion(1,ia)+Gy(n)*Rion(2,ia)+Gz(n)*Rion(3,ia)
          rbox16=rbox16+conjg(rhoe_G(n))*dVloc_G(n,Kion(ia))*exp(-zI*Gd)
        end do
      end do
      Ebox1(3)=Ebox1(3)+rbox13/dble(NPUW)
      Ebox1(4)=Ebox1(4)+rbox14/dble(NPUW)
      Ebox1(5)=Ebox1(5)+rbox15/dble(NPUW)
      Ebox1(6)=Ebox1(6)+rbox16/dble(NPUW)
    end do
  end do
end select

sum_temp3(1:4)=Ebox1(3:6)
call comm_summation(sum_temp3,sum_temp4,4,nproc_group_global)
Ebox2(3:6)=sum_temp4(1:4)


elp3(1405)=get_wtime()
elp3(1454)=elp3(1454)+elp3(1405)-elp3(1404)

do iik=k_sta,k_end
do iob=1,iobnum
!$OMP parallel
!$OMP do private(iatom,jj)
  do iatom=1,MI
  do jj=1,maxlm
    uVpsibox3(jj,iatom,iob,iik)=0.d0
    uVpsibox4(jj,iatom,iob,iik)=0.d0
  end do
  end do
!$OMP end parallel
end do
end do

! 7:Enl_l
do iatom=1,MI
  ik=Kion(iatom)
  do jj=1,Mps(iatom)
    x=(dble(Jxyz(1,jj,iatom)-1)-dble(Jxxyyzz(1,jj,iatom)*lg_num(1)))*Hgs(1)
    y=(dble(Jxyz(2,jj,iatom)-1)-dble(Jxxyyzz(2,jj,iatom)*lg_num(2)))*Hgs(2)
    z=(dble(Jxyz(3,jj,iatom)-1)-dble(Jxxyyzz(3,jj,iatom)*lg_num(3)))*Hgs(3)
    do iik=k_sta,k_end
      ekr(jj,iatom,iik)=exp(zi*(k_rd(1,iik)*x+k_rd(2,iik)*y+k_rd(3,iik)*z))
    end do
  end do
end do

do iik=k_sta,k_end
do iob=1,iobnum
  do iatom=1,MI
    ik=Kion(iatom)
    loop_lm3 : do lm=1,(Mlps(ik)+1)**2
      if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm3
      sumbox=0.d0
        do jj=1,Mps(iatom)
          sumbox=sumbox+uV(jj,lm,iatom)*tzpsi_in(Jxyz(1,jj,iatom),Jxyz(2,jj,iatom),Jxyz(3,jj,iatom),iob,iik)*ekr(jj,iatom,iik)
        end do
      uVpsibox3(lm,iatom,iob,iik)=sumbox*Hvol
    end do loop_lm3
  end do
end do
end do


call comm_summation(uVpsibox3,uVpsibox4,maxlm*MI*iobnum*k_num,nproc_group_korbital)

elp3(1406)=get_wtime()
elp3(1455)=elp3(1455)+elp3(1406)-elp3(1405)

do iik=k_sta,k_end
do iob=1,iobnum
  call calc_allob(iob,p_allob)
  do iatom=1,MI
  ik=Kion(iatom)
  loop_lm4 : do lm=1,(Mlps(ik)+1)**2
    if ( abs(uVu(lm,iatom))<1.d-5 ) cycle loop_lm4
    Ebox1(7)=Ebox1(7)+rocc(p_allob,iik)*wtk(iik)*abs(uVpsibox4(lm,iatom,iob,iik))**2/uVu(lm,iatom)
  end do loop_lm4
  end do
end do
end do

sum_temp1=Ebox1(7)
call comm_summation(sum_temp1,sum_temp2,nproc_group_rho)
Ebox2(7)=sum_temp2

elp3(1407)=get_wtime()
elp3(1456)=elp3(1456)+elp3(1407)-elp3(1406)

Etot=sum(Ebox2(1:7))+Exc

end if

return

end subroutine Total_Energy_periodic_scf
