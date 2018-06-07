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
!-----------------------------------------------------------------------------------------
module builtin_pz_sp
  implicit none
  
contains
  
  
  subroutine exc_cor_pz_sp(nl, rho_s, exc, eexc, vexc_s)
    implicit none
    integer, intent(in) :: nl
    real(8), intent(in) :: rho_s(nl,2)
    real(8), intent(out) :: exc(nl), eexc(nl), vexc_s(nl,2)
    integer :: i
    real(8) :: trho, trho_s(2), e_xc, vxc_s(2)

    do i=1,NL
      trho=rho_s(i,1)+rho_s(i,2)+1.d-20
      trho_s(1:2)=2*rho_s(i,1:2)+5.d-21
      call PZxc_sp(trho,trho_s,e_xc,vxc_s)
      exc(i)=e_xc
      Eexc(i)=e_xc*trho
      vexc_s(i,1:2)=vxc_s(1:2)
    enddo
    return
  end subroutine exc_cor_pz_sp



  subroutine PZxc_sp(trho,trho_s,exc,vxc_s)
    implicit none
    real(8),parameter :: Pi=3.141592653589793d0
    real(8),parameter :: gammaU=-0.1423d0,gammaP=-0.0843d0
    real(8),parameter :: beta1U=1.0529d0,beta2U=0.3334d0
    real(8),parameter :: beta1P=1.3981d0,beta2P=0.2611d0
    real(8),parameter :: AU=0.0311d0,BU=-0.048d0
    real(8),parameter :: CU=0.002d0,DU=-0.0116d0
    real(8),parameter :: AP=0.01555d0,BP=-0.0269d0
    real(8),parameter :: CP=0.0007d0,DP=-0.0048d0
    real(8),parameter :: const = 0.75d0*(3d0/(2d0*pi))**(2d0/3d0)
    real(8) :: rs,cx
    real(8) :: zeta,sf,dsf
    real(8) :: ecu,ecp,vcu,vcp
    real(8) :: sgnsigma(2)
    integer :: js
    real(8),intent(in) :: trho,trho_s(2)
    real(8),intent(out) :: exc
    real(8),intent(out) :: vxc_s(2)

    sgnsigma(1)=1.d0
    sgnsigma(2)=-1.d0
    rs=(3.d0/(4.d0*Pi*trho))**(1.d0/3.d0)

    zeta=(trho_s(1)-trho_s(2))/trho
    sf=((1.d0+zeta)**(4.d0/3.d0)+(1.d0-zeta)**(4.d0/3.d0)-2.d0)/(2.d0**(4.d0/3.d0)-2.d0)
! dsf=d(smallf)/d(zeta)
    dsf=4.d0/3.d0*((1.d0+zeta)**(1.d0/3.d0)-(1.d0-zeta)**(1.d0/3.d0))/(2.d0**(4.d0/3.d0)-2.d0)

    cx=3.d0/4.d0*(3.d0/(2*Pi))**(2.d0/3.d0)

    exc=-cx/rs+(-2.d0**(1.d0/3.d0)+1.d0)*cx/rs*sf
    do js=1,2
      vxc_s(js)=exc-1.d0/3.d0*cx/rs+1.d0/3.d0*(-2.d0**(1.d0/3.d0)+1.d0)*cx/rs*sf &
                             +(-2.d0**(1.d0/3.d0)+1.d0)*cx/rs*(sgnsigma(js)-zeta)*dsf
    end do

    if ( rs > 1.d0 ) then
      ecu=gammaU/(1.d0+beta1U*sqrt(rs)+beta2U*rs)
      ecp=gammaP/(1.d0+beta1P*sqrt(rs)+beta2P*rs)
      vcu=gammaU*(1.d0+7.d0/6.d0*beta1U*sqrt(rs)+4.d0/3.d0*beta2U*rs)/(1.d0+beta1U*sqrt(rs)+beta2U*rs)**2
      vcp=gammaP*(1.d0+7.d0/6.d0*beta1P*sqrt(rs)+4.d0/3.d0*beta2P*rs)/(1.d0+beta1P*sqrt(rs)+beta2P*rs)**2
    else 
      Ecu=AU*log(rs)+BU+CU*rs*log(rs)+DU*rs
      Ecp=AP*log(rs)+BP+CP*rs*log(rs)+DP*rs
      Vcu=AU*log(rs)+(BU-AU/3.d0)+2.d0/3.d0*CU*rs*log(rs)+(2.d0*DU-CU)/3.d0*rs
      Vcp=AP*log(rs)+(BP-AP/3.d0)+2.d0/3.d0*CP*rs*log(rs)+(2.d0*DP-CP)/3.d0*rs
    end if
  
    exc=exc+ecu+sf*(ecp-ecu)
    do js=1,2
      vxc_s(js)=vxc_s(js)+vcu+sf*(vcp-vcu)+(ecp-ecu)*(sgnsigma(js)-zeta)*dsf
    end do

    return
  end Subroutine PZxc_sp
  
end module builtin_pz_sp
