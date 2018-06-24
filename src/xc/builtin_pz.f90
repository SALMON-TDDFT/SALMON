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
module builtin_pz
  implicit none
  
contains
  
  
  subroutine exc_cor_pz(nl, rho_s, exc, eexc, vexc)
    implicit none
    integer, intent(in) :: nl
    real(8), intent(in) :: rho_s(nl)
    real(8), intent(out) :: exc(nl), eexc(nl), vexc(nl)
    integer :: i
    real(8) :: trho, e_xc, de_xc_drho

    ! call rho_j_tau(gs_rt, rho_s, tau_s, j_s, grho_s, lrho_s)
    ! rho_s=rho*0.5d0
    ! if(flag_nlcc)rho_s = rho_s + 0.5d0*rho_nlcc
    
    !!$omp parallel do private(i, trho, rs, v_xc, e_xc, rssq, rsln)
    do i=1,NL
      trho=2*rho_s(i)
      call PZxc(trho,e_xc,de_xc_drho)
      exc(i)=e_xc
      Eexc(i)=e_xc*trho
      Vexc(i)=e_xc+trho*de_xc_drho
    enddo
    return
  end subroutine exc_cor_pz



  Subroutine PZxc(trho,exc,dexc_drho)
    implicit none
    real(8),parameter :: Pi=3.141592653589793d0
    real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
    real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
    real(8),parameter :: CU=0.002d0,DU=-0.0116d0
    real(8),parameter :: const = 0.75d0*(3d0/(2d0*pi))**(2d0/3d0)
    real(8) :: ttrho,rs,rssq,rsln
    real(8),intent(in) :: trho
    real(8),intent(out) :: exc
    real(8),intent(out) :: dexc_drho


    ttrho=trho+1d-10
    rs=(3d0/(4*Pi*ttrho))**(1d0/3d0)
    exc=-const/rs
    dexc_drho=exc/(3*ttrho)
    if (rs>1d0) then
      rssq=sqrt(rs)
      exc=exc+gammaU/(1d0+beta1U*rssq+beta2U*rs)
      dexc_drho=dexc_drho+gammaU*(0.5d0*beta1U*rssq+beta2U*rs)/(3*ttrho)/(1+beta1U*rssq+beta2U*rs)**2
    else
      rsln=log(rs)
      exc=exc+AU*rsln+BU+CU*rs*rsln+DU*rs
      dexc_drho=dexc_drho-rs/(3*ttrho)*(AU/rs+CU*(rsln+1)+DU)
    endif
    return
  End Subroutine PZxc
  
end module builtin_pz
