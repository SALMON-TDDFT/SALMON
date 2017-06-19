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
MODULE Total_Energy_sub

use scf_data
use hpsi2_sub
use new_world_sub
use read_pslfile_sub
use inner_product_sub

INTERFACE Total_Energy
  MODULE PROCEDURE R_Total_Energy, C_Total_Energy
END INTERFACE

CONTAINS
!=======================================================================

SUBROUTINE R_Total_Energy(psi_in)
use salmon_parallel, only: nproc_group_global, nproc_group_orbital, nproc_group_h
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
implicit none

real(8) :: psi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),  &
                1:iobnum,1)
real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd)
real(8) :: htpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
integer :: iob,ia,ib
integer :: ix,iy,iz
real(8) :: rab
real(8) :: esp2(itotMST,1)
real(8) :: sum1,sum2
real(8) :: rbox
integer :: iob_allob

iwk_size=2
call make_iwksta_iwkend

ihpsieff=0

esp=0.d0

!$OMP parallel do
do iz=mg_sta(3)-Nd,mg_end(3)+Nd
do iy=mg_sta(2)-Nd,mg_end(2)+Nd
do ix=mg_sta(1)-Nd,mg_end(1)+Nd
  tpsi(ix,iy,iz)=0.d0
end do
end do
end do

do iob=1,iobnum
  call calc_allob(iob,iob_allob)

!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    tpsi(ix,iy,iz)=psi_in(ix,iy,iz,iob,1)
  end do
  end do
  end do

  call hpsi2(tpsi,htpsi,iob_allob,0,0)

  rbox=0.d0
!$OMP parallel do reduction ( + : rbox )
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rbox=rbox+psi_in(ix,iy,iz,iob,1)*htpsi(ix,iy,iz)
  end do
  end do
  end do

  esp(iob_allob,1)=dble(rbox)*Hvol

end do

call comm_summation(esp,esp2,itotMST,nproc_group_global)
esp=esp2

Etot=0.d0
Etot = Etot + sum( rocc(:itotMST,1)*esp(:itotMST,1) )*wtk(1)

do ia=1,MI
do ib=1,ia-1
  rab=sqrt((Rion(1,ia)-Rion(1,ib))**2      &
           +(Rion(2,ia)-Rion(2,ib))**2      &
           +(Rion(3,ia)-Rion(3,ib))**2)
  Etot=Etot+Zps(Kion(ia))*Zps(Kion(ib))/rab
end do
end do

if(ilsda == 0)then
  sum1=0.d0
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    sum1=sum1+rho(ix,iy,iz)*(-0.5d0*Vh(ix,iy,iz)-Vxc(ix,iy,iz))
  end do
  end do
  end do
  call comm_summation(sum1,sum2,nproc_group_h)
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
  call comm_summation(sum1,sum2,nproc_group_orbital)
  Etot=Etot+sum2*Hvol+Exc
end if

return

END SUBROUTINE R_Total_Energy

!=======================================================================

SUBROUTINE C_Total_Energy(psi_in)
use salmon_parallel, only: nproc_group_global, nproc_group_orbital, nproc_group_h
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
implicit none

complex(8) :: psi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),  &
                1:iobnum,1)
complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd)
complex(8) :: htpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
integer :: iob,ia,ib
integer :: ix,iy,iz
real(8) :: rab
real(8) :: esp2(itotMST,1)
real(8) :: sum1,sum2
complex(8) :: cbox
integer :: iob_allob

elp3(861)=get_wtime()

iwk_size=2
call make_iwksta_iwkend

ihpsieff=0

esp2=0.d0

!$OMP parallel do
do iz=mg_sta(3)-Nd,mg_end(3)+Nd
do iy=mg_sta(2)-Nd,mg_end(2)+Nd
do ix=mg_sta(1)-Nd,mg_end(1)+Nd
  tpsi(ix,iy,iz)=0.d0
end do
end do
end do

elp3(862)=get_wtime()
elp3(882)=elp3(882)+elp3(862)-elp3(861)
do iob=1,iobnum
  call calc_allob(iob,iob_allob)
  elp3(863)=get_wtime()

!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    tpsi(ix,iy,iz)=psi_in(ix,iy,iz,iob,1)
  end do
  end do
  end do

  call hpsi2(tpsi,htpsi,iob_allob,0,0)

  elp3(864)=get_wtime()
  elp3(884)=elp3(884)+elp3(864)-elp3(863)

  cbox=0.d0
!$OMP parallel do reduction ( + : cbox )
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    cbox=cbox+conjg(psi_in(ix,iy,iz,iob,1))*htpsi(ix,iy,iz)
  end do
  end do
  end do
  
  esp2(iob_allob,1)=dble(cbox)*Hvol

  elp3(865)=get_wtime()
  elp3(885)=elp3(885)+elp3(865)-elp3(864)
end do
elp3(866)=get_wtime()

call comm_summation(esp2,esp,itotMST,nproc_group_global)

elp3(867)=get_wtime()
elp3(887)=elp3(887)+elp3(867)-elp3(866)

Etot=0.d0
Etot = Etot + sum( rocc(:itotMST,1)*esp(:itotMST,1) )*wtk(1)

elp3(868)=get_wtime()
elp3(888)=elp3(888)+elp3(868)-elp3(867)

do ia=1,MI
do ib=1,ia-1
  rab=sqrt((Rion(1,ia)-Rion(1,ib))**2      &
           +(Rion(2,ia)-Rion(2,ib))**2      &
           +(Rion(3,ia)-Rion(3,ib))**2)
  Etot=Etot+Zps(Kion(ia))*Zps(Kion(ib))/rab
end do
end do

elp3(869)=get_wtime()
elp3(889)=elp3(889)+elp3(869)-elp3(868)

if(ilsda == 0)then
  sum1=0.d0
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    sum1=sum1+rho(ix,iy,iz)*(-0.5d0*Vh(ix,iy,iz)-Vxc(ix,iy,iz))
  end do
  end do
  end do
  call comm_summation(sum1,sum2,nproc_group_h)
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
  call comm_summation(sum1,sum2,nproc_group_orbital)
  Etot=Etot+sum2*Hvol+Exc
end if

elp3(870)=get_wtime()
elp3(890)=elp3(890)+elp3(870)-elp3(869)
elp3(891)=elp3(891)+elp3(870)-elp3(861)

return

END SUBROUTINE C_Total_Energy

END MODULE Total_Energy_sub
