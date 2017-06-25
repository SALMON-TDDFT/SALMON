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
!=================================================== LDA (Perdew-Zunger)

SUBROUTINE Exc_Cor_ns
use scf_data

implicit none
integer :: ix,iy,iz
real(8) :: trho(mg_sta(1):mg_end(1),  &
                mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3))
real(8) :: trho_s(mg_sta(1):mg_end(1),  &
                mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3),2)

!$OMP parallel do private(iz,iy,ix) collapse(3)
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  if(rho(ix,iy,iz)>=0.d0)then
    trho(ix,iy,iz)=rho(ix,iy,iz)+1.d-20
  else if(rho(ix,iy,iz)<-1.d-5)then
    write(*,*) "rho has large negative value", ix,iy,iz,rho(ix,iy,iz)
    stop
  else
    trho(ix,iy,iz)=1.d-20
  end if
end do
end do
end do
if(ilsda == 1) then
!$OMP parallel do private(iz,iy,ix) collapse(3)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    trho_s(ix,iy,iz,1:2)=rho_s(ix,iy,iz,1:2)+5.d-21
    trho(ix,iy,iz)=rho(ix,iy,iz)+1.d-20
  end do
  end do
  end do
end if

call xc_LDA(trho,trho_s)

END SUBROUTINE Exc_Cor_ns

!======================================================================

SUBROUTINE xc_LDA(trho,trho_s)
use salmon_parallel, only: nproc_group_h
use salmon_communication, only: comm_summation
use scf_data
use new_world_sub
implicit none
integer :: js

real(8),parameter :: gm=-0.1423d0
real(8),parameter :: gmu=-0.1423d0,gmp=-0.0843d0
! letter s means small, letter l means large
real(8),parameter :: sbu1=1.0529d0, sbu2=0.3334d0
real(8),parameter :: sbp1=1.3981d0, sbp2=0.2611d0
real(8),parameter :: rlau=0.0311d0 , rlbu=-0.048d0
real(8),parameter :: rlcu=0.002d0 , rldu=-0.0116d0
real(8),parameter :: rlap=0.01555d0 , rlbp=-0.0269d0
real(8),parameter :: rlcp=0.0007d0 , rldp=-0.0048d0

real(8) :: sgnsigma(2)
real(8) :: Cx
real(8) :: rs
real(8) :: zeta
real(8) :: sf
real(8) :: dsf

real(8),allocatable :: Ecu(:,:,:),Ecp(:,:,:)
real(8),allocatable :: Vcu(:,:,:),Vcp(:,:,:)

integer :: ix,iy,iz
real(8) :: sum1

real(8) :: trho(mg_sta(1):mg_end(1),  &
                mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3))
real(8) :: trho_s(mg_sta(1):mg_end(1),  &
                mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3),2)

real(8) :: Ex(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3))
real(8) :: Ec(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3))
real(8) :: Ex_LDA(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3))
real(8) :: Ec_LDA(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
             mg_sta(3):mg_end(3))

real(8) :: Vx(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
               mg_sta(3):mg_end(3))
real(8) :: Vc(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
               mg_sta(3):mg_end(3))
real(8) :: Vx_LDA(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                   mg_sta(3):mg_end(3))
real(8) :: Vc_LDA(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                   mg_sta(3):mg_end(3))
real(8) :: Vxc_LDA(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                    mg_sta(3):mg_end(3))
real(8) :: Vx_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                 mg_sta(3):mg_end(3),2)
real(8) :: Vc_s(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                 mg_sta(3):mg_end(3),2)
real(8) :: Vx_s_LDA(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                     mg_sta(3):mg_end(3),2)
real(8) :: Vc_s_LDA(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                     mg_sta(3):mg_end(3),2)
real(8) :: Vxc_s_LDA(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
                      mg_sta(3):mg_end(3),2)

if(ilsda==1)then
  allocate(Ecu(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
               mg_sta(3):mg_end(3)))
  allocate(Ecp(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
               mg_sta(3):mg_end(3)))
  allocate(Vcu(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
               mg_sta(3):mg_end(3)))
  allocate(Vcp(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),  &
               mg_sta(3):mg_end(3)))
end if

sgnsigma(1)=1.d0
sgnsigma(2)=-1.d0

!!$OMP parallel do &
!!$OMP private(Cx,rs,zeta,sf,dsf,iz,iy,ix) collapse(3)
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  
  rs=(3.d0/(4.d0*Pi*trho(ix,iy,iz)))**(1.d0/3.d0)
  
  if(ilsda == 1) then
     zeta=(trho_s(ix,iy,iz,1)-trho_s(ix,iy,iz,2))/trho(ix,iy,iz)
  
     sf=((1.d0+zeta)**(4.d0/3.d0)+(1.d0-zeta)**(4.d0/3.d0)-2.d0)      &
                                         /(2.d0**(4.d0/3.d0)-2.d0)
! dsf=d(smallf)/d(zeta)
     dsf=4.d0/3.d0*((1.d0+zeta)**(1.d0/3.d0) &
                   -(1.d0-zeta)**(1.d0/3.d0)) &
                  /(2.d0**(4.d0/3.d0)-2.d0)
  end if
  
! Exchange
  
  Cx=3.d0/4.d0*(3.d0/(2*Pi))**(2.d0/3.d0)
  
  if(ilsda == 0)then
    Ex_LDA(ix,iy,iz)=-Cx/rs
    Vx_LDA(ix,iy,iz)=-4.d0/3.d0*Cx/rs
  else if(ilsda == 1)then
    Ex_LDA(ix,iy,iz)=-Cx/rs+(-2.d0**(1.d0/3.d0)+1.d0)*Cx/rs*sf
    do js=1,2
      Vx_s_LDA(ix,iy,iz,js)=Ex_LDA(ix,iy,iz)-1.d0/3.d0*Cx/rs  &
                             +1.d0/3.d0*(-2.d0**(1.d0/3.d0)+1.d0)      &
                                  *Cx/rs*sf &
                             +(-2.d0**(1.d0/3.d0)+1.d0)*Cx/rs &
                                  *(sgnsigma(js)-zeta)*dsf
    end do
  end if
  
  
! Correlation
! Ceperley and Alder exchange correlation written in J. P. Perdew and A. Zunger, Phys. Rev. B, vol. 23, 5048 (1981).
  
  if(ilsda == 0)then
    if ( rs > 1.d0 ) then
      Ec_LDA(ix,iy,iz)=gm/( 1.d0 + sbu1*sqrt(rs) + sbu2*rs )
      Vc_LDA(ix,iy,iz)=gm*( 1.d0 + 7.d0/6.d0*sbu1*sqrt(rs) + 4.d0/3.d0*sbu2*rs ) &
                     /( 1.d0 + sbu1*sqrt(rs) + sbu2*rs )**2
    else 
      Ec_LDA(ix,iy,iz)=rlau*log(rs) + rlbu   &
                     + rlcu*rs*log(rs) + rldu*rs
      Vc_LDA(ix,iy,iz)=rlau*log(rs) + (rlbu-rlau/3.d0)   &
                     + 2.d0/3.d0*rlcu*rs*log(rs)      &
                     + (2.d0*rldu-rlcu)/3.d0*rs
    end if
  else if(ilsda == 1)then
    if ( rs > 1.d0 ) then
      Ecu(ix,iy,iz)=gmu/( 1.d0 + sbu1*sqrt(rs) + sbu2*rs )
      Ecp(ix,iy,iz)=gmp/( 1.d0 + sbp1*sqrt(rs) + sbp2*rs )

      Vcu(ix,iy,iz)=gmu*( 1.d0 + 7.d0/6.d0*sbu1*sqrt(rs) + 4.d0/3.d0*sbu2*rs )      &
            /( 1.d0 + sbu1*sqrt(rs) + sbu2*rs )**2
      Vcp(ix,iy,iz)=gmp*( 1.d0 + 7.d0/6.d0*sbp1*sqrt(rs) + 4.d0/3.d0*sbp2*rs )      &
            /( 1.d0 + sbp1*sqrt(rs) + sbp2*rs )**2
    else 
      Ecu(ix,iy,iz)=rlau*log(rs) + rlbu + rlcu*rs*log(rs) + rldu*rs
      Ecp(ix,iy,iz)=rlap*log(rs) + rlbp + rlcp*rs*log(rs) + rldp*rs

      Vcu(ix,iy,iz)=rlau*log(rs) + (rlbu-rlau/3.d0) + 2.d0/3.d0*rlcu*rs*log(rs)      &
                      + (2.d0*rldu-rlcu)/3.d0*rs
      Vcp(ix,iy,iz)=rlap*log(rs) + (rlbp-rlap/3.d0) + 2.d0/3.d0*rlcp*rs*log(rs)      &
                      + (2.d0*rldp-rlcp)/3.d0*rs
    end if
  
    Ec_LDA(ix,iy,iz)=Ecu(ix,iy,iz)+sf*(Ecp(ix,iy,iz)-Ecu(ix,iy,iz))
    do js=1,2
      Vc_s_LDA(ix,iy,iz,js)=Vcu(ix,iy,iz) +sf*(Vcp(ix,iy,iz)-Vcu(ix,iy,iz))      &
                                +(Ecp(ix,iy,iz)-Ecu(ix,iy,iz))      &
                                 *(sgnsigma(js)-zeta)*dsf
    end do
  end if
  
  if(ilsda==0)then
    Vxc_LDA(ix,iy,iz)=(Vx_LDA(ix,iy,iz)+Vc_LDA(ix,iy,iz))
  else if(ilsda==1)then
    Vxc_s_LDA(ix,iy,iz,1:2)=(Vx_s_LDA(ix,iy,iz,1:2)+Vc_s_LDA(ix,iy,iz,1:2))
  end if
  
  if(ilsda==0)then
    Ex(ix,iy,iz)=Ex_LDA(ix,iy,iz)
    Ec(ix,iy,iz)=Ec_LDA(ix,iy,iz)
    Vx(ix,iy,iz)=Vx_LDA(ix,iy,iz)
    Vc(ix,iy,iz)=Vc_LDA(ix,iy,iz)
    Vxc(ix,iy,iz)=Vxc_LDA(ix,iy,iz)
  else if(ilsda==1)then
    Ex(ix,iy,iz)=Ex_LDA(ix,iy,iz)
    Ec(ix,iy,iz)=Ec_LDA(ix,iy,iz)
    Vx_s(ix,iy,iz,1:2)=Vx_s_LDA(ix,iy,iz,1:2)
    Vc_s(ix,iy,iz,1:2)=Vc_s_LDA(ix,iy,iz,1:2)
    Vxc_s(ix,iy,iz,1:2)=Vxc_s_LDA(ix,iy,iz,1:2)
  end if
  
end do
end do
end do

sum1=0.d0
!$omp parallel do reduction(+ : sum1) private(iz,iy,ix) collapse(3)
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  sum1=sum1+trho(ix,iy,iz)*(Ex(ix,iy,iz)+Ec(ix,iy,iz))
end do
end do
end do
sum1=sum1*Hvol
call comm_summation(sum1,Exc,nproc_group_h)

if(ilsda==1) deallocate(Ecu,Ecp,Vcu,Vcp)

END SUBROUTINE xc_LDA

!======================================================================
