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
SUBROUTINE Total_Energy_groupob(tzpsi_in,htpsi,ifunc)
use salmon_parallel, only: nproc_group_global, nproc_group_h, nproc_group_korbital
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use scf_data
use new_world_sub
use inner_product_sub
use read_pslfile_sub
implicit none

complex(8) :: tzpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd,  &
                1:iobnum,k_sta:k_end)
complex(8) :: htpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd,  &
                1:iobnum,k_sta:k_end)
complex(8) :: tzpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd,  &
                1:iobnum,k_sta:k_end)
integer :: ifunc
integer :: iob,ia,ib,iik
integer :: ix,iy,iz
real(8) :: rab
real(8) :: sum1,sum2
real(8) :: rbox
complex(8) :: cbox

if(ifunc==1)then

  iwk_size=2
  call make_iwksta_iwkend
  
  ihpsieff=0
  
  esp2=0.d0
  
  if(ilsda==1)then 
!$OMP parallel do private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Vlocal(ix,iy,iz,1)=Vpsl(ix,iy,iz)+Vh(ix,iy,iz)+Vxc_s(ix,iy,iz,1)
      Vlocal(ix,iy,iz,2)=Vpsl(ix,iy,iz)+Vh(ix,iy,iz)+Vxc_s(ix,iy,iz,2)
    end do
    end do
    end do
  end if
  
  
  call hpsi_groupob(tzpsi_in,htpsi,tzpsi_out,Vlocal,0,0)
  
  
  do iik=k_sta,k_end
  do iob=1,iobnum
    cbox=0.d0
!$OMP parallel do reduction ( + : cbox ) private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      cbox=cbox+conjg(tzpsi_in(ix,iy,iz,iob,iik))*htpsi(ix,iy,iz,iob,iik)
    end do
    end do
    end do
  
    esp2(iob,iik)=dble(cbox)*Hvol
  
  end do
  end do
  
  elp3(761)=get_wtime()
  call comm_summation(esp2,esp,itotMST*num_kpoints_rd,nproc_group_global)
  
  
  elp3(762)=get_wtime()
  elp3(782)=elp3(782)+elp3(762)-elp3(761)

else if(ifunc==2)then

  elp3(761)=get_wtime()
  call comm_summation(esp2,esp,itotMST*num_kpoints_rd,nproc_group_global)

  elp3(762)=get_wtime()
  elp3(782)=elp3(782)+elp3(762)-elp3(761)

end if


if(iflag_md==0)then
  Etot=Eion
else if(iflag_md==1)then
  Etot=0.d0
  do ia=1,MI
  do ib=1,ia-1
    rab=sqrt((Rion(1,ia)-Rion(1,ib))**2      &
            +(Rion(2,ia)-Rion(2,ib))**2      &
            +(Rion(3,ia)-Rion(3,ib))**2)
    Etot=Etot+Zps(Kion(ia))*Zps(Kion(ib))/rab
  end do
  end do
end if

do iik=k_sta,k_end
  rbox=0.d0
!$OMP parallel do reduction ( + : rbox )
  do iob=1,itotMST
    rbox = rbox + rocc(iob,iik)*esp(iob,iik) *wtk(iik)
  end do
  Etot=Etot+rbox
end do

if(ilsda == 0)then
  if((ifunc==1.and.mod(itt,2)==1).or.(ifunc==2.and.mod(itt,2)==0))then
    sum1=0.d0
!$OMP parallel do reduction (+ : sum1 ) private(iz,iy,ix)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      sum1=sum1+rho(ix,iy,iz)*(-0.5d0*Vh_stock2(ix,iy,iz)-Vxc(ix,iy,iz))
    end do
    end do
    end do
  else
    sum1=0.d0
!$OMP parallel do reduction (+ : sum1 ) private(iz,iy,ix)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      sum1=sum1+rho(ix,iy,iz)*(-0.5d0*Vh_stock1(ix,iy,iz)-Vxc(ix,iy,iz))
    end do
    end do
    end do
  end if
  elp3(761)=get_wtime()
  call comm_summation(sum1,sum2,nproc_group_h)
  elp3(762)=get_wtime()
  elp3(783)=elp3(783)+elp3(762)-elp3(761)
  Etot=Etot+sum2*Hvol+Exc
else if(ilsda == 1)then
  sum1=0.d0
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    sum1=sum1+rho(ix,iy,iz)*(-0.5d0*Vh(ix,iy,iz))   &
             -Vxc_s(ix,iy,iz,1)*rho_s(ix,iy,iz,1)   &
             -Vxc_s(ix,iy,iz,2)*rho_s(ix,iy,iz,2) 
  end do
  end do
  end do
  call comm_summation(sum1,sum2,nproc_group_korbital)
  Etot=Etot+sum2*Hvol+Exc
end if

return

END SUBROUTINE Total_Energy_groupob

