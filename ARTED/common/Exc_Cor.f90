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
!This file is "Exc_Cor.f90"
subroutine Exc_Cor(GS_RT,NBtmp,zu)
  use Global_Variables
  use timer
  implicit none
  integer,intent(in)       :: GS_RT
  integer,intent(in)       :: NBtmp
  complex(8),intent(inout) :: zu(NL,NBtmp,NK_s:NK_e)
  call timer_begin(LOG_EXC_COR)
  if(functional == 'PZ')    call Exc_Cor_PZ
  if(functional == 'PZM')   call Exc_Cor_PZM
  if(functional == 'PBE')   call Exc_Cor_PBE(GS_RT)
  if(functional == 'TBmBJ') call Exc_Cor_TBmBJ(GS_RT)
  if(functional == 'TPSS')  call Exc_Cor_TPSS(GS_RT)
  if(functional == 'VS98')  call Exc_Cor_VS98(GS_RT)
  call timer_end(LOG_EXC_COR)

contains
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Exc_Cor_PZ()
  use Global_Variables
  implicit none
  real(8) :: rho_s(NL)
  integer :: i
  real(8) :: trho,e_xc,de_xc_drho

!$acc kernels
!  call rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
  rho_s=rho*0.5d0
  if(flag_nlcc)rho_s = rho_s + 0.5d0*rho_nlcc
!!$omp parallel do private(i,trho,rs,V_xc,E_xc,rssq,rsln)
  do i=1,NL
    trho=2*rho_s(i)
    call PZxc(trho,e_xc,de_xc_drho)
    Eexc(i)=e_xc*trho
    Vexc(i)=e_xc+trho*de_xc_drho
  enddo
!$acc end kernels
  return
End Subroutine Exc_Cor_PZ
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Exc_Cor_PZM()
  use Global_Variables
  implicit none
  real(8) :: rho_s(NL)
  integer :: i
  real(8) :: trho,e_xc,de_xc_drho

!$acc kernels
!  call rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
  rho_s=rho*0.5d0
  if(flag_nlcc)rho_s = rho_s + 0.5d0*rho_nlcc
!!$omp parallel do private(i,trho,rs,V_xc,E_xc,rssq,rsln)
  do i=1,NL
    trho=2*rho_s(i)
    call PZMxc(trho,e_xc,de_xc_drho)
    Eexc(i)=e_xc*trho
    Vexc(i)=e_xc+trho*de_xc_drho
  enddo
!$acc end kernels
  return
End Subroutine Exc_Cor_PZM
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Exc_Cor_PBE(GS_RT)
  use Global_Variables
  implicit none
  integer,intent(in) :: GS_RT
  real(8) :: rho_s(NL),tau_s(NL),j_s(NL,3),grho_s(NL,3),lrho_s(NL)
  real(8) :: agrho_s(NL)
  integer :: i
  real(8) :: omega
  real(8) :: FxHSE,drho_FxHSE,ds_FxHSE,ec_unif,drs_ec_unif,H,drs_H,dt_H
  real(8) :: rho_t,rs,kf,ks,grho(NL,3),agrho(NL),srho,trho
  real(8) :: exlda,Ex(NL),Ec(NL),Vx(NL),Vc(NL)
  real(8) :: dFxdrho_tbl(NL,3),dHksrho_tbl(NL,3)

  call rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
  agrho_s(:)=sqrt(grho_s(:,1)**2+grho_s(:,2)**2+grho_s(:,3)**2)
  grho(:,1)=2*grho_s(:,1)
  grho(:,2)=2*grho_s(:,2)
  grho(:,3)=2*grho_s(:,3)
  agrho(:)=2*agrho_s(:)
  omega=0d0

  do i=1,NL
    rho_t=rho(i)
    if(rho(i) < 1d-10) rho_t=rho(i)+1d-10
    rs=(3d0/(4*pi*rho_t))**(1d0/3d0)
    kf=(3*pi**2*rho_t)**(1d0/3d0)
    ks=sqrt(4*kf/pi)
    srho=agrho(i)/(2*kf*rho_t)
    trho=agrho(i)/(2*ks*rho_t)
    call HSEFx(rho_t,srho,omega,FxHSE,drho_FxHSE,ds_FxHSE)
    call PBEc_old(rho_t,trho,ec_unif,drs_ec_unif,H,drs_H,dt_H)
    exlda=-3d0/(4*pi)*kf
    Ex(i)=exlda*FxHSE
    Vx(i)=exlda*(4d0/3d0+rho(i)*drho_FxHSE-4*srho/3d0*ds_FxHSE)
    Ec(i)=ec_unif+H
    Vc(i)=ec_unif+H-rs/3*(drs_ec_unif+drs_H)-7d0/6d0*trho*dt_H
    dFxdrho_tbl(i,:)=ds_FxHSE*grho(i,:)/agrho(i)
    dHksrho_tbl(i,:)=dt_H/(2*ks)*grho(i,:)/agrho(i)
  enddo

  do i=1,NL
    Vx(i)=Vx(i)+3d0/(8*pi)*( &
 &    nabx(1)*(dFxdrho_tbl(ifdx(1,i),1)-dFxdrho_tbl(ifdx(-1,i),1)) &
 &   +nabx(2)*(dFxdrho_tbl(ifdx(2,i),1)-dFxdrho_tbl(ifdx(-2,i),1)) &
 &   +nabx(3)*(dFxdrho_tbl(ifdx(3,i),1)-dFxdrho_tbl(ifdx(-3,i),1)) &
 &   +nabx(4)*(dFxdrho_tbl(ifdx(4,i),1)-dFxdrho_tbl(ifdx(-4,i),1)) & 
 &   +naby(1)*(dFxdrho_tbl(ifdy(1,i),2)-dFxdrho_tbl(ifdy(-1,i),2)) &
 &   +naby(2)*(dFxdrho_tbl(ifdy(2,i),2)-dFxdrho_tbl(ifdy(-2,i),2)) &
 &   +naby(3)*(dFxdrho_tbl(ifdy(3,i),2)-dFxdrho_tbl(ifdy(-3,i),2)) &
 &   +naby(4)*(dFxdrho_tbl(ifdy(4,i),2)-dFxdrho_tbl(ifdy(-4,i),2)) & 
 &   +nabz(1)*(dFxdrho_tbl(ifdz(1,i),3)-dFxdrho_tbl(ifdz(-1,i),3)) &
 &   +nabz(2)*(dFxdrho_tbl(ifdz(2,i),3)-dFxdrho_tbl(ifdz(-2,i),3)) &
 &   +nabz(3)*(dFxdrho_tbl(ifdz(3,i),3)-dFxdrho_tbl(ifdz(-3,i),3)) &
 &   +nabz(4)*(dFxdrho_tbl(ifdz(4,i),3)-dFxdrho_tbl(ifdz(-4,i),3)) )
    Vc(i)=Vc(i)-( &
 &    nabx(1)*(dHksrho_tbl(ifdx(1,i),1)-dHksrho_tbl(ifdx(-1,i),1)) &
 &   +nabx(2)*(dHksrho_tbl(ifdx(2,i),1)-dHksrho_tbl(ifdx(-2,i),1)) &
 &   +nabx(3)*(dHksrho_tbl(ifdx(3,i),1)-dHksrho_tbl(ifdx(-3,i),1)) &
 &   +nabx(4)*(dHksrho_tbl(ifdx(4,i),1)-dHksrho_tbl(ifdx(-4,i),1)) & 
 &   +naby(1)*(dHksrho_tbl(ifdy(1,i),2)-dHksrho_tbl(ifdy(-1,i),2)) &
 &   +naby(2)*(dHksrho_tbl(ifdy(2,i),2)-dHksrho_tbl(ifdy(-2,i),2)) &
 &   +naby(3)*(dHksrho_tbl(ifdy(3,i),2)-dHksrho_tbl(ifdy(-3,i),2)) &
 &   +naby(4)*(dHksrho_tbl(ifdy(4,i),2)-dHksrho_tbl(ifdy(-4,i),2)) & 
 &   +nabz(1)*(dHksrho_tbl(ifdz(1,i),3)-dHksrho_tbl(ifdz(-1,i),3)) &
 &   +nabz(2)*(dHksrho_tbl(ifdz(2,i),3)-dHksrho_tbl(ifdz(-2,i),3)) &
 &   +nabz(3)*(dHksrho_tbl(ifdz(3,i),3)-dHksrho_tbl(ifdz(-3,i),3)) &
 &   +nabz(4)*(dHksrho_tbl(ifdz(4,i),3)-dHksrho_tbl(ifdz(-4,i),3)) )
  enddo
  Eexc=rho*(Ex+Ec)
  Vexc=Vx+Vc
  return
End Subroutine Exc_Cor_PBE
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Exc_Cor_TBmBJ(GS_RT)
  use Global_Variables
  implicit none
  real(8),parameter :: alpha=-0.012d0,beta=1.023d0,gamma=0.80d0
  integer,intent(in) :: GS_RT
  real(8) :: rho_s(NL),tau_s(NL),j_s(NL,3),grho_s(NL,3),lrho_s(NL)
  real(8) :: c,tau_s_jrho,D_s_jrho,Q_s,rhs,x_s,b_s,Vx_BR,Vx_MBJ
  real(8) :: trho,rs,rhos,ec,dec_drhoa,dec_drhob
  integer :: i

  call rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
  if(cval > 1d-10) then
    c=cval ! use c-value given by input file
  else
    c=sum(sqrt(grho_s(:,1)**2+grho_s(:,2)**2+grho_s(:,3)**2)/rho_s(:))*Hxyz/aLxyz
    c=alpha+beta*sqrt(c)
  endif

  do i=1,NL
    tau_s_jrho=tau_s(i)-(j_s(i,1)**2+j_s(i,2)**2+j_s(i,3)**2)/rho_s(i)/2
    D_s_jrho=2*tau_s_jrho-0.25d0*(grho_s(i,1)**2+grho_s(i,2)**2+grho_s(i,3)**2)/rho_s(i)
    Q_s=(lrho_s(i)-2*gamma*D_s_jrho)/6d0
    rhs=2d0/3d0*pi**(2d0/3d0)*rho_s(i)**(5d0/3d0)/Q_s
    call BR_Newton(rhs,x_s)
    b_s=(x_s**3*exp(-x_s)/(8*pi*rho_s(i)))**(1d0/3d0)
    Vx_BR=-(1-exp(-x_s)-0.5d0*x_s*exp(-x_s))/b_s
    Vx_MBJ=c*Vx_BR+(3*c-2)/pi*sqrt(5d0/12d0)*sqrt(2*tau_s_jrho/rho_s(i))

    trho=rho(i)+1d-10
    rs=(3d0/(4*Pi*trho))**(1d0/3d0)
    Eexc(i)=-.4582d0/rs
    Vexc(i)=Vx_MBJ
    rhos=rho_s(i)
    call PWc(rhos,rhos,ec,dec_drhoa,dec_drhob)
    Vexc(i)=Vexc(i)+ec+rho(i)*dec_drhoa
    Eexc(i)=Eexc(i)+ec
    Eexc(i)=Eexc(i)*trho
  enddo
  return
End Subroutine Exc_Cor_TBmBJ
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
SUBROUTINE BR_Newton(rhs,x_s)
  implicit none
  real(8),intent(IN) :: rhs
  real(8),intent(OUT) :: x_s
  integer iter
  real(8) :: xmin,xmax,x,fx,dfx

! find xmax
  xmin=0d0
  x=1d0
  do
    fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
    if(fx < 0) then
      xmax=x
      exit
    endif
    x=x*2
  enddo
! bi-section
  do
    x=0.5d0*(xmin+xmax)
    fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
    if(fx < 0) then
      xmax=x
    else
      xmin=x
    endif
    if(xmax-xmin < 1d-4) exit
  enddo
! Newton-Raphson
  do iter=1,5
    fx=x*exp(-2d0/3d0*x)/rhs-(x-2)
    dfx=(1-2d0/3d0*x)*exp(-2d0/3d0*x)/rhs-1d0
    x=x-fx/dfx
  enddo
  x_s=x
  return
end subroutine BR_Newton
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Exc_Cor_TPSS(GS_RT)
  use Global_Variables
  implicit none
  integer,intent(in) :: GS_RT
  integer i
  real(8),parameter :: kappa=0.804d0,c=1.59096d0,e=1.537d0,mu=0.21951d0
  real(8),parameter :: b=0.40d0,d=2.8d0,C0=0.53d0
  real(8) :: rho_s(NL),tau_s(NL),j_s(NL,3),grho_s(NL,3),lrho_s(NL)
  real(8) :: agrho_s(NL),e_xc(NL),dexc_drho(NL),dexc_dtau(NL),dexc_dgrho(NL,3)
  real(8) :: trho,agrho,tau,grho(3),p,z,alp,qb,sqe,x,Fx,ex_unif,ex,dx_dqb
  real(8) :: dalp_dp,dqb_dalp,dqb_dp,dx_dp,dFx_dp,dalp_dz,dqb_dz,dx_dz,dFx_dz
  real(8) :: dz_dtaua,ec_up,dec_drhoa_up,dec_drhob_up,dec_dgrhoa_up(3)
  real(8) :: dec_dgrhob_up(3),zero,gzero(3),ec_p,dec_drhoa_p,dec_drhob_p
  real(8) :: dec_dgrhoa_p(3),dec_dgrhob_p(3),ec_tilde,dec_tilde_drhoa
  real(8) :: dec_tilde_dgrhoa(3),ec_rev,dec_rev_drhoa,dec_rev_dgrhoa(3)
  real(8) :: dec_rev_dtaua,dex_drho,dex_dgrho(3),dex_dtau
  real(8) :: ec,dec_drho,dec_dgrho(3),dec_dtau

  call rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
  agrho_s(:)=sqrt(grho_s(:,1)**2+grho_s(:,2)**2+grho_s(:,3)**2)
  tau_s(:)=tau_s(:)-(j_s(:,1)**2+j_s(:,2)**2+j_s(:,3)**2)/(2*rho_s(:))

  do i=1,NL
    trho=rho_s(i)*2
    agrho=agrho_s(i)*2
    tau=tau_s(i)*2
    grho(:)=grho_s(i,:)*2
    p=agrho**2/(4*(3*Pi**2)**(2d0/3d0)*trho**(8d0/3d0))
    z=agrho**2/(8*trho)/tau
    alp=5*p/3*(1/z-1)
    qb=9d0/20d0*(alp-1)/(1+b*alp*(alp-1))**0.5d0+2d0*p/3d0
    sqe=sqrt(e)
    x=((10d0/81d0+c*z**2/(1+z**2)**2)*p+146d0/2025d0*qb**2 &
 &     -73d0/405d0*qb*sqrt((0.6d0*z)**2/2d0+p**2/2)+(10d0/81d0)**2*p**2/kappa &
 &     +2*sqe*10d0/81d0*(0.6d0*z)**2+e*mu*p**3)/(1+sqe*p)**2
    Fx=1+kappa-kappa/(1+x/kappa)
    ex_unif=-3d0/(4d0*Pi)*(3*Pi**2*trho)**(1d0/3d0)
    ex=trho*ex_unif*Fx

    dx_dqb=(146d0/2025d0*2*qb-73d0/405d0*sqrt((0.6d0*z)**2/2+p**2/2))/(1+sqe*p)**2
    dalp_dp=5d0/3d0*(1/z-1)
    dqb_dalp=9d0/20d0/(1+b*alp*(alp-1))**0.5d0-0.5d0*9d0/20d0*(alp-1)*b*(2*alp-1)/(1+b*alp*(alp-1))**1.5d0
    dqb_dp=dalp_dp*dqb_dalp+2d0/3d0
    dx_dp=((10d0/81d0+c*z**2/(1+z**2)**2)-73d0/405d0*qb*(p/2)/sqrt((0.6d0*z)**2/2+p**2/2) &
 &        +(10d0/81d0)**2/kappa*2*p+3*e*mu*p**2)/(1+sqe*p)**2 &
 &       -2*sqe*x/(1+sqrt(e)*p)+dx_dqb*dqb_dp
    dFx_dp=dx_dp/(1+x/kappa)**2
    dalp_dz=-5*p/(3*z**2)
    dqb_dz=dalp_dz*dqb_dalp
    dx_dz=(2*c*p*z*(1-z**2)/(1+z**2)**3+2*sqrt(e)*10d0/81d0*0.36d0*2*z)/(1+sqe*p)**2 &
 &       +dx_dqb*dqb_dz
    dFx_dz=dx_dz/(1+x/kappa)**2
    dex_drho=4d0/3d0*ex_unif*Fx-ex_unif*(8d0/3d0*p*dFx_dp+z*dFx_dz)

    dex_dgrho(:)=trho*ex_unif*(dFx_dp*(2*grho(:))/(4*(3*Pi**2)**(2d0/3d0)*trho**(8d0/3d0)) &
 &                            +dFx_dz*grho(:)/(4*trho*tau))
    dz_dtaua=-z/tau
    dex_dtau=trho*ex_unif*dz_dtaua*dFx_dz

    call PBEc(rho_s(i),rho_s(i),grho_s(i,:),grho_s(i,:),ec_up,dec_drhoa_up,dec_drhob_up,dec_dgrhoa_up,dec_dgrhob_up)
    zero=0d0
    gzero(:)=0d0
    call PBEc(rho_s(i),zero    ,grho_s(i,:),gzero(:),   ec_p, dec_drhoa_p ,dec_drhob_p, dec_dgrhoa_p, dec_dgrhob_p )

    ec_tilde=ec_up
    dec_tilde_drhoa=dec_drhoa_up
    dec_tilde_dgrhoa(:)=dec_dgrhoa_up(:)
    if(ec_p > ec_up) then
      ec_tilde=ec_p
      dec_tilde_drhoa=dec_drhoa_p/2
      dec_tilde_dgrhoa(:)=dec_dgrhoa_p(:)/2
    endif
    ec_rev=ec_up*(1+C0*z**2)-(1+C0)*z**2*ec_tilde
    ec=trho*ec_rev*(1+d*ec_rev*z**3)

    dec_rev_drhoa=(1+C0*z**2)*dec_drhoa_up-2*C0*z**2/trho*ec_up+(1+C0)*2*z**2/trho*ec_tilde-(1+C0)*z**2*dec_tilde_drhoa
    dec_drho=ec_rev-2*d*z**3*ec_rev**2+trho*dec_rev_drhoa*(1+2*d*z**3*ec_rev)

    dec_rev_dgrhoa(:)=(1+C0*z**2)*dec_dgrhoa_up(:)+2*C0*z*grho(:)/(4*trho*tau)*ec_up &
 &     -(1+C0)*2*z*grho(:)/(4*trho*tau)*ec_tilde-(1+C0)*z**2*dec_tilde_dgrhoa(:)
    dec_dgrho(:)=trho*dec_rev_dgrhoa(:)*(1+2*d*z**3*ec_rev)+grho(:)*3*d*z**2/(4*tau)*ec_rev**2

    dec_rev_dtaua=-2*z**2/tau*(C0*ec_up-(1+C0)*ec_tilde)
    dec_dtau=trho*dec_rev_dtaua*(1+2*d*z**3*ec_rev)-3*trho*d*z**3/tau*ec_rev**2

    e_xc(i)=ex+ec
    dexc_drho(i)=dex_drho+dec_drho
    dexc_dtau(i)=dex_dtau+dec_dtau
    dexc_dgrho(i,:)=dex_dgrho(:)+dec_dgrho(:)
  enddo

  tmass(:)=dexc_dtau(:)
  tjr(:,1)=tmass(:)*j_s(:,1)/rho_s(:)
  tjr(:,2)=tmass(:)*j_s(:,2)/rho_s(:)
  tjr(:,3)=tmass(:)*j_s(:,3)/rho_s(:)
  tjr2(:)=tmass*(j_s(:,1)**2+j_s(:,2)**2+j_s(:,3)**2)/rho_s(:)**2

  Eexc(:)=rho(:)*e_xc(:)
  do i=1,NL
    Vexc(i)=dexc_drho(i)-( &
 &    nabx(1)*(dexc_dgrho(ifdx(1,i),1)-dexc_dgrho(ifdx(-1,i),1)) &
 &   +nabx(2)*(dexc_dgrho(ifdx(2,i),1)-dexc_dgrho(ifdx(-2,i),1)) &
 &   +nabx(3)*(dexc_dgrho(ifdx(3,i),1)-dexc_dgrho(ifdx(-3,i),1)) &
 &   +nabx(4)*(dexc_dgrho(ifdx(4,i),1)-dexc_dgrho(ifdx(-4,i),1)) & 
 &   +naby(1)*(dexc_dgrho(ifdy(1,i),2)-dexc_dgrho(ifdy(-1,i),2)) &
 &   +naby(2)*(dexc_dgrho(ifdy(2,i),2)-dexc_dgrho(ifdy(-2,i),2)) &
 &   +naby(3)*(dexc_dgrho(ifdy(3,i),2)-dexc_dgrho(ifdy(-3,i),2)) &
 &   +naby(4)*(dexc_dgrho(ifdy(4,i),2)-dexc_dgrho(ifdy(-4,i),2)) & 
 &   +nabz(1)*(dexc_dgrho(ifdz(1,i),3)-dexc_dgrho(ifdz(-1,i),3)) &
 &   +nabz(2)*(dexc_dgrho(ifdz(2,i),3)-dexc_dgrho(ifdz(-2,i),3)) &
 &   +nabz(3)*(dexc_dgrho(ifdz(3,i),3)-dexc_dgrho(ifdz(-3,i),3)) &
 &   +nabz(4)*(dexc_dgrho(ifdz(4,i),3)-dexc_dgrho(ifdz(-4,i),3)) )
  enddo

  return
  end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Exc_Cor_VS98(GS_RT)
  use Global_Variables
  implicit none
  integer,intent(in) :: GS_RT
  real(8),parameter :: CF=3d0/5d0*(3*Pi**2)**(2d0/3d0)
  real(8),parameter :: alp_ex=0.00186726d0,alp_aa=0.00515088d0,alp_ab=0.00304966d0
  real(8),parameter :: a_ex=-0.9800683d0,a_aa=0.3270912d0,a_ab=0.7035010d0
  real(8),parameter :: b_ex=-0.003556788d0,b_aa=-0.03228915d0,b_ab=0.007694574d0
  real(8),parameter :: c_ex=0.006250326d0,c_aa=-0.02942406d0,c_ab=0.05152765d0
  real(8),parameter :: d_ex=-0.00002354518d0,d_aa=0.002134222d0,d_ab=0.00003394308d0
  real(8),parameter :: e_ex=-0.0001282732d0,e_aa=-0.005451559d0,e_ab=-0.001269420d0
  real(8),parameter :: f_ex=0.0003574822d0,f_aa=0.01577575d0,f_ab=0.001296118d0
  integer :: i
  real(8) :: x_s,z_s,x,z,D_s,zero,fxz_ex,dfxz_ex_dx,dfxz_ex_dz
  real(8) :: fxz_ab,dfxz_ab_dx,dfxz_ab_dz,fxz_aa,dfxz_aa_dx,dfxz_aa_dz
  real(8) :: rhos,rhos13,rhos43,rhos53,ec,dec_drhoa,dec_drhob
  real(8) :: ecp,decp_drhoa,decp_drhob,ec_ab,dec_ab_drhoa,ec_aa,dec_aa_drhoa
  real(8) :: rho_s(NL),tau_s(NL),j_s(NL,3),grho_s(NL,3),lrho_s(NL)
  real(8) :: agrho_s(NL),dex_drho(NL),dex_dtau(NL),dex_dgrho(NL,3)
  real(8) :: dec_drho(NL),dec_dtau(NL),dec_dgrho(NL,3),dexc_dgrho(NL,3)

  call rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
  agrho_s(:)=sqrt(grho_s(:,1)**2+grho_s(:,2)**2+grho_s(:,3)**2)
  tau_s(:)=tau_s(:)-(j_s(:,1)**2+j_s(:,2)**2+j_s(:,3)**2)/(2*rho_s(:))

  do i=1,NL
    rhos=rho_s(i)
    rhos13=rhos**(1d0/3d0)
    rhos43=rhos**(4d0/3d0)
    rhos53=rhos**(5d0/3d0)
    x_s=agrho_s(i)/rhos43
    z_s=tau_s(i)/rhos53-CF
    call fec_xz(x_s,z_s,alp_ex,a_ex,b_ex,c_ex,d_ex,e_ex,f_ex,fxz_ex,dfxz_ex_dx,dfxz_ex_dz)
    Eexc(i)=2*rhos43*fxz_ex
    dex_drho(i)=rhos13*(4d0/3d0*fxz_ex-4d0/3d0*x_s*dfxz_ex_dx &
 &                                   -5d0/3d0*(z_s+CF)*dfxz_ex_dz)
    dex_dtau(i)=dfxz_ex_dz/rhos13
    dex_dgrho(i,:)=grho_s(i,:)/agrho_s(i)*dfxz_ex_dx

    x=sqrt(2*x_s**2)
    z=z_s*2
    D_s=1-x_s**2/(4*(z_s+CF))
    call fec_xz(x,z,alp_ab,a_ab,b_ab,c_ab,d_ab,e_ab,f_ab,fxz_ab,dfxz_ab_dx,dfxz_ab_dz)
    call fec_xz(x_s,z_s,alp_aa,a_aa,b_aa,c_aa,d_aa,e_aa,f_aa,fxz_aa,dfxz_aa_dx,dfxz_aa_dz)
    zero=0d0
    call PWc(rhos,rhos,ec,dec_drhoa,dec_drhob)
    call PWc(rhos,zero,ecp,decp_drhoa,decp_drhob)
    ec_ab=rhos*2*(ec-ecp)
    dec_ab_drhoa=ec+rhos*2*dec_drhoa-ecp-rhos*decp_drhoa
    ec_aa=rhos*ecp
    dec_aa_drhoa=ecp+rhos*decp_drhoa
    Eexc(i)=Eexc(i)+fxz_ab*ec_ab+fxz_aa*D_s*ec_aa*2
    dec_drho(i)=(-4d0/3d0*x_s**2/(x*rhos)*dfxz_ab_dx-5/(3d0*rhos)*(z_s+CF)*dfxz_ab_dz)*ec_ab &
 &             +fxz_ab*dec_ab_drhoa+(-4d0/3d0*x_s/rhos*dfxz_aa_dx-5d0/(3d0*rhos)*(z_s+CF)*dfxz_aa_dz)*D_s*ec_aa &
 &             +fxz_aa/rhos*(1-D_s)*ec_aa+fxz_aa*D_s*dec_aa_drhoa
    dec_dtau(i)=dfxz_ab_dz/rhos53*ec_ab+dfxz_aa_dz/rhos53*D_s*ec_aa &
 &             +fxz_aa*x_s**2/4d0/rhos53/(z_s+CF)**2*ec_aa
    dec_dgrho(i,:)=grho_s(i,:)/agrho_s(i)*(x_s/x/rhos43*dfxz_ab_dx*ec_ab &
 &                                        +(dfxz_aa_dx*D_s-x_s/2d0/(z_s+CF)*fxz_aa)/rhos43*ec_aa)
  enddo
  dexc_dgrho=dex_dgrho+dec_dgrho
  tmass=dex_dtau+dec_dtau
  tjr(:,1)=tmass(:)*j_s(:,1)/rho_s(:)
  tjr(:,2)=tmass(:)*j_s(:,2)/rho_s(:)
  tjr(:,3)=tmass(:)*j_s(:,3)/rho_s(:)
  tjr2(:)=tmass*(j_s(:,1)**2+j_s(:,2)**2+j_s(:,3)**2)/rho_s(:)**2

  do i=1,NL
    Vexc(i)=dex_drho(i)+dec_drho(i)-( &
 &    nabx(1)*(dexc_dgrho(ifdx(1,i),1)-dexc_dgrho(ifdx(-1,i),1)) &
 &   +nabx(2)*(dexc_dgrho(ifdx(2,i),1)-dexc_dgrho(ifdx(-2,i),1)) &
 &   +nabx(3)*(dexc_dgrho(ifdx(3,i),1)-dexc_dgrho(ifdx(-3,i),1)) &
 &   +nabx(4)*(dexc_dgrho(ifdx(4,i),1)-dexc_dgrho(ifdx(-4,i),1)) & 
 &   +naby(1)*(dexc_dgrho(ifdy(1,i),2)-dexc_dgrho(ifdy(-1,i),2)) &
 &   +naby(2)*(dexc_dgrho(ifdy(2,i),2)-dexc_dgrho(ifdy(-2,i),2)) &
 &   +naby(3)*(dexc_dgrho(ifdy(3,i),2)-dexc_dgrho(ifdy(-3,i),2)) &
 &   +naby(4)*(dexc_dgrho(ifdy(4,i),2)-dexc_dgrho(ifdy(-4,i),2)) & 
 &   +nabz(1)*(dexc_dgrho(ifdz(1,i),3)-dexc_dgrho(ifdz(-1,i),3)) &
 &   +nabz(2)*(dexc_dgrho(ifdz(2,i),3)-dexc_dgrho(ifdz(-2,i),3)) &
 &   +nabz(3)*(dexc_dgrho(ifdz(3,i),3)-dexc_dgrho(ifdz(-3,i),3)) &
 &   +nabz(4)*(dexc_dgrho(ifdz(4,i),3)-dexc_dgrho(ifdz(-4,i),3)) )
  enddo
  return
  end subroutine

SUBROUTINE fec_xz(x,z,alp,a,b,c,d,e,f,fxz,dfxz_dx,dfxz_dz)
  implicit none
  real(8),intent(IN) :: x,z,alp,a,b,c,d,e,f
  real(8),intent(OUT) :: fxz,dfxz_dx,dfxz_dz
  real(8) :: gamma

  gamma=1+alp*(x**2+z)
  fxz=a/gamma+(b*x**2+c*z)/gamma**2+(d*x**4+e*x**2*z+f*z**2)/gamma**3
  dfxz_dx=2*alp*x*(-a/gamma**2-2*(b*x**2+c*z)/gamma**3-3*(d*x**4+e*x**2*z+f*z**2)/gamma**4) &
 &      +2*b*x/gamma**2+(4*d*x**3+2*e*x*z)/gamma**3
  dfxz_dz=alp*(-a/gamma**2-2*(b*x**2+c*z)/gamma**3-3*(d*x**4+e*x**2*z+f*z**2)/gamma**4) &
 &      +c/gamma**2+(e*x**2+2*f*z)/gamma**3
  return
  end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine PZxc(trho,exc,dexc_drho)
!$acc routine seq
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine PZMxc(trho,exc,dexc_drho)
!$acc routine seq
  implicit none
  real(8),parameter :: Pi=3.141592653589793d0
  real(8),parameter :: gammaU=-0.1423d0,beta1U=1.0529d0
  real(8),parameter :: beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
  real(8),parameter :: CU=0.2019151940622859d-2, DU=-0.1163206637891297d-1
  real(8),parameter :: const = 0.75d0*(3d0/(2d0*pi))**(2d0/3d0)
  real(8) :: ttrho,rs,rssq,rsln
  real(8),intent(in) :: trho
  real(8),intent(out) :: exc
  real(8),intent(out) :: dexc_drho

  ttrho=trho+1d-10
  rs=(3d0/(4*Pi*ttrho))**(1d0/3d0)
  exc=-const/rs
  dexc_drho=exc/(3d0*ttrho)
  if (rs>1d0) then
    rssq=sqrt(rs)
    exc=exc+gammaU/(1d0+beta1U*rssq+beta2U*rs)
    dexc_drho=dexc_drho+gammaU*(0.5d0*beta1U*rssq+beta2U*rs)/(3d0*ttrho)/(1d0+beta1U*rssq+beta2U*rs)**2
  else
    rsln=log(rs)
    exc=exc+AU*rsln+BU+CU*rs*rsln+DU*rs
    dexc_drho=dexc_drho-rs/(3d0*ttrho)*(AU/rs+CU*(rsln+1d0)+DU)
  endif
  return
End Subroutine PZMxc
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
  SUBROUTINE HSEFx(rho,ss,omega,FxHSE,drho_FxHSE,ds_FxHSE)
  implicit none

  real(8),intent(IN) :: rho,ss,omega
  real(8),intent(OUT) :: FxHSE,drho_FxHSE,ds_FxHSE

  real(8),parameter :: pi=3.1415926535897932d0
  real(8),parameter :: sqpi=1.7724538509055160d0
  real(8),parameter :: eps=1d-10
  real(8),parameter :: A=1.0161144D0,B=-3.7170836D-1
  real(8),parameter :: C=-7.7215461D-2,D=5.7786348D-1
  real(8),parameter :: E=-5.1955731D-2
  real(8),parameter :: Ha1=9.79681d-3,Ha2=4.10834d-2
  real(8),parameter :: Ha3=1.87440d-1,Ha4=1.20824d-3
  real(8),parameter :: Ha5=3.47188d-2
  real(8),parameter :: Fc1=6.4753871d0,Fc2=4.7965830d-1
  real(8),parameter :: ea1=-1.128223946706117d0,ea2=1.452736265762971d0
  real(8),parameter :: ea3=-1.243162299390327d0,ea4=0.971824836115601d0
  real(8),parameter :: ea5=-0.568861079687373d0,ea6=0.246880514820192d0
  real(8),parameter :: ea7=-0.065032363850763d0,ea8=0.008401793031216d0
  real(8),parameter :: wcut=14d0
  real(8),parameter :: EGscut=0.08d0,EGa1=-2.628417880d-2
  real(8),parameter :: EGa2=-7.117647788d-2,EGa3=8.534541323d-2
  real(8),parameter :: expcut=700d0
  real(8),parameter :: exei1=4.03640D0,exei2=1.15198D0
  real(8),parameter :: exei3=5.03627D0,exei4=4.19160D0
  real(8),parameter :: smax=8.572844D0,strans=8.3d0
  real(8),parameter :: sconst=1.879622316D1

  real(8) s,s2,s3,s4,s5,s6,kf,w,w2,w3,w4,w5,w6,w7,w8
  real(8) H,ds_H,ds_Hs,F,ds_CFs,eb1,Hsbw,Hsbw2,Hsbw3,Hsbw4
  real(8) DHsb,DHsb2,DHsb3,DHsb4,DHsb5,DHsb12,DHsb32,DHsb52,DHsb72,DHsb92
  real(8) A94,A942,A943,A944,HA94,HA942,HA943,HA945
  real(8) DHs,DHs2,DHs3,DHs4,DHs52,DHs72
  real(8) DHsw,DHsw12,DHsw32,DHsw52,DHsw72
  real(8) wsqDHsw,wsqDHsw3,wsqDHsw5,wsqDHsw7
  real(8) G_a,G_b,EG,ds_EGs,exer,exei
!  real(8) ei
!  real(8) derfc

  s=ss
! s is modified to enforce Lieb-Oxford bound
  if(s.gt.strans) then
    s=smax-sconst/s**2
  endif

  s2=s**2; s3=s**3; s4=s**4; s5=s**5; s6=s**6

  kf=(3d0*pi**2*rho)**(1d0/3d0)
  w=omega/kf

  w2=w**2; w3=w**3; w4=w**4; w5=w**5; w6=w**6; w7=w**7; w8=w**8

  H=(Ha1*s2+Ha2*s4)/(1d0+Ha3*s4+Ha4*s5+Ha5*s6)
  ds_H=(2*Ha1*s+4*Ha2*s3)/(1+Ha3*s4+Ha4*s5+Ha5*s6) &
 &    -(Ha1*s2+Ha2*s4)/(1+Ha3*s4+Ha4*s5+Ha5*s6)**2 &
 &    *(4*Ha3*s3+5*Ha4*s4+6*Ha5*s5)
  ds_Hs=ds_H*s2+2*s*H

  F=Fc1*H+Fc2
  ds_CFs=C*(2*s*F+s2*Fc1*ds_H)

  if(w.lt.wcut) then
    eb1=1.455915450052607d0
  else
    eb1=2.0d0
  endif

  Hsbw=s2*H+eb1*w2
  Hsbw2=Hsbw**2; Hsbw3=Hsbw**3; Hsbw4=Hsbw**4

  DHsb=D+s2*H+eb1*w2
  DHsb2=DHsb**2; DHsb3=DHsb**3; DHsb4=DHsb**4; DHsb5=DHsb**5
  DHsb12=DHsb**0.5d0; DHsb32=DHsb**1.5d0; DHsb52=DHsb**2.5d0
  DHsb72=DHsb**3.5d0; DHsb92=DHsb**4.5d0

  A94=9d0/(4d0*A)
  A942=A94**2; A943=A94**3; A944=A94**4

  HA94=A94*Hsbw
  HA942=HA94**2; HA943=HA94**3; HA945=HA94**5

  DHs=D+s2*H
  DHs2=DHs**2; DHs3=DHs**3; DHs4=DHs**4
  DHs52=DHs**2.5d0; DHs72=DHs**3.5d0

  DHsw=DHs+w2
  DHsw12=DHsw**0.5d0; DHsw32=DHsw**1.5d0; DHsw52=DHsw**2.5d0; DHsw72=DHsw**3.5d0

  wsqDHsw=w/DHsw12
  wsqDHsw3=wsqDHsw**3; wsqDHsw5=wsqDHsw**5; wsqDHsw7=wsqDHsw**7

  if(s.gt.EGscut) then
    G_a=sqpi*(15d0*E+6d0*C*(1d0+F*s2)*DHs+4d0*B*DHs2+8d0*A*DHs3) &
 &     /(16d0*DHs72) &
 &     -0.75d0*pi*sqrt(A)*exp(A94*H*s2) &
 &     *derfc(1.5d0*s*sqrt(H/A))
    G_b=(15d0/16d0*sqpi*s2)/DHs72
    EG=-(0.75d0*pi+G_a)/G_b
    ds_EGs=-14d0/5d0*sqpi*DHs52*ds_Hs &
 &         +sqpi*sqrt(A)/5*DHs52  &
 &         *exp(9*H*s2/(4*A))*derfc(1.5d0*s*sqrt(H/A)) &
 &         *(14+9/A*DHs)*ds_Hs &
 &         -6d0/5d0*DHs72/sqrt(H*s2)*ds_Hs &
 &         -2d0/5d0*DHs*ds_CFs-2d0/5d0*C*(1+F*s2)*ds_Hs &
 &         -8d0/15d0*B*DHs*ds_Hs-8d0/5d0*A*DHs2*ds_Hs
  else
    EG=EGa1+EGa2*s2+EGa3*s4
    ds_EGs=2*s*EG+s2*(2*EGa2*s+4*EGa3*s3)
  endif

  if(HA94 > eps) then
    if(HA94 < expcut) then
      exer=pi*exp(HA94)*derfc(sqrt(HA94))
      exei=exp(HA94)*ei(-HA94)
    else
      exer=pi*(1d0/(sqpi*sqrt(HA94))-1d0/(2d0*sqrt(pi*HA943)) &
 &    +3d0/(4d0*sqrt(pi*HA945)))
      exei=-(1d0/HA94) &
 &        *(HA942+exei1*HA94+exei2)/(HA942+exei3*HA94+exei4)
    endif
  endif

!  if(w < eps.and.s < eps) then
  if(HA94 < eps) then
    FxHSE=1d0
    drho_FxHSE=0d0
    ds_FxHSE=0d0
    return
  endif

  if(w < eps) then
    FxHSE=-8d0/9d0*(-A/2*(log(DHs)-log(H*s2))-A/2*exei &
 &                  +B/2/DHs+C*(1+F*s2)/2/DHs**2+(E+EG*s2)/DHs3)
    drho_FxHSE=0d0
    ds_FxHSE=-8d0/9d0*(ds_CFs/2/DHs2+ds_EGs/DHs3 &
 &                    +ds_Hs*(A/2*(-1/DHs-A94*exei)-B/2/DHs2-C*(1+F*s2)/DHs3-3*(E+EG*s2)/DHs4))
    return
  endif

  FxHSE=-8d0/9d0*( &
 & -A/2*(log(DHsb)-log(Hsbw))-A/2*exei &
 & +B/2/DHs*(1-wsqDHsw) &
 & +3d0/4d0*C*(1+F*s2)/DHs2*(2d0/3d0-wsqDHsw+wsqDHsw3/3) &
 & +15d0/8d0*(E+EG*s2)/DHs3*(8d0/15d0-wsqDHsw+2d0/3d0*wsqDHsw3-wsqDHsw5/5) &
 & )

  drho_FxHSE= -8d0/9d0*(2*w/3d0/sqpi/rho)*( &
 &   -3d0/4d0*sqrt(A)*pi*exp(A94*(H*s2+w2))*derfc(sqrt(A94*(H*s2+w2))) &
 &   +sqpi*(A/2/DHsw12+B/4/DHsw32+C*(1+F*s2)*3d0/8d0/DHsw52+(E+EG*s2)*15d0/16d0/DHsw72) )

  ds_FxHSE= -8d0/9d0*( &
 &  3d0/4d0*ds_CFs/DHs2*(2d0/3d0-wsqDHsw+wsqDHsw3/3) &
 & +15d0/8d0*ds_EGs/DHs3*(8d0/15d0-wsqDHsw+2d0/3d0*wsqDHsw3-wsqDHsw5/5) &
 & -ds_Hs*( &
 &         9d0/8d0*exei+A/2/DHsb &
 &        +3d0/4d0*B/DHs2*(2d0/3d0-wsqDHsw+wsqDHsw3/3) &
 &        +15d0/8d0*C*(1+F*s2)/DHs3*(8d0/15d0-wsqDHsw+2d0/3d0*wsqDHsw3-wsqDHsw5/5) &
 &        +105d0/16d0*(E+EG*s2)/DHs4*(16d0/35d0-wsqDHsw+wsqDHsw3 &
 &                                    -3d0/5d0*wsqDHsw5+wsqDHsw7/7) &
 &          ))

  if(w > wcut) return

  FxHSE=FxHSE-8d0/9d0*( &
 &   -9d0/8d0*exei*(-ea2*w2+A94*ea4*w4-A942*ea6*w6+A943*ea8*w8) &
 &   -9d0/8d0*(ea4*w4/Hsbw+ea6*w6*(-A94/Hsbw+1/Hsbw2)+ea8*w8*(A942/Hsbw-A94/Hsbw2+2/Hsbw3)) &
 &   -9d0/8d0*exer/sqrt(A94)*(ea1*w-ea3*w3*A94+ea5*w5*A942-ea7*w7*A943) &
 &   -9d0/8d0*sqpi/sqrt(Hsbw)*(ea3*w3+ea5*w5*(-A94+1/(2*Hsbw)) &
 &                 +ea7*w7*(A942-A94/Hsbw/2+3/(4*Hsbw2))) &
 &   +A/2*(ea2*w2/DHsb+ea4*w4/DHsb2+ea6*w6*2/DHsb3+ea8*w8*6/DHsb4) &
 &   +A*sqpi*(ea1*w/2/DHsb12+ea3*w3/4/DHsb32+ea5*w5*3d0/8d0/DHsb52+ea7*w7*15d0/16d0/DHsb72) )

  ds_FxHSE=ds_FxHSE-8d0/9d0*( &
 &   -ds_Hs*(-9d0/8d0*((ea2*w2*A94-ea4*w4*A942+ea6*w6*A943-ea8*w8*A944)*exei &
 &                     +ea2*w2/Hsbw+ea4*w4*(-A94/Hsbw+1/Hsbw2) &
 &                     +ea6*w6*(A942/Hsbw-A94/Hsbw2+2/Hsbw3) &
 &                     +ea8*w8*(-A943/Hsbw+A942/Hsbw2-2*A94/Hsbw3+6/Hsbw4)) &
 &           -9d0/8d0*(sqrt(A94)*(-ea1*w+ea3*w3*A94-ea5*w5*A942+ea7*w7*A943)*exer &
 &                    +sqpi/sqrt(Hsbw)*(ea1*w+ea3*w3*(-A94+1/Hsbw/2) &
 &                                     +ea5*w5*(A942-A94/2/Hsbw+3d0/4d0/Hsbw2) &
 &                                     +ea7*w7*(-A943+A942/2/Hsbw &
 &                                     -3d0/4d0*A94/Hsbw2+15d0/8d0/Hsbw3))) &
 &           +A/2*(ea2*w2/DHsb2+ea4*w4*2/DHsb3+ea6*w6*6/DHsb4+ea8*w8*24/DHsb5) &
 &           +A*sqpi*(ea1*w/4/DHsb32+ea3*w3*3d0/8d0/DHsb52 &
 &                   +ea5*w5*15d0/16d0/DHsb72+ea7*w7*105d0/32d0/DHsb92) ))
  end subroutine

  FUNCTION ei(x)
  implicit none
  real(8) ei,x,xx
  integer,parameter :: one=1
  if(x.ge.0d0) stop 'bad x in ei'
  xx=-x
  ei=-expint(one,xx)
  return
  end function

  FUNCTION expint(n,x)
  implicit none
  INTEGER n,MAXIT
  REAL(8) expint,x,EPS,FPMIN,EULER
  PARAMETER (MAXIT=100,EPS=1.d-7,FPMIN=1.d-30,EULER=.5772156649d0)
  INTEGER i,ii,nm1
  REAL(8) a,b,c,d,del,fact,h,psi
  nm1=n-1
  if(n.lt.0.or.x.lt.0.d0.or.(x.eq.0.d0.and.(n.eq.0.or.n.eq.1))) then
    stop 'bad arguments in expint'
  else if(n.eq.0)then
    expint=exp(-x)/x
  else if(x.eq.0.d0)then
    expint=1.d0/nm1
  else if(x.gt.1.d0)then
    b=x+n
    c=1.d0/FPMIN
    d=1.d0/b
    h=d
    do i=1,MAXIT
      a=-i*(nm1+i)
      b=b+2.d0
      d=1.d0/(a*d+b)
      c=b+a/c
      del=c*d
      h=h*del
      if(abs(del-1.d0).lt.EPS)then
        expint=h*exp(-x)
        return
      endif
    enddo
    stop 'continued fraction failed in expint'
  else
    if(nm1.ne.0)then
      expint=1.d0/nm1
    else
      expint=-log(x)-EULER
    endif
    fact=1.d0
    do i=1,MAXIT
      fact=-fact*x/i
      if(i.ne.nm1) then
        del=-fact/(i-nm1)
      else
        psi=-EULER
        do ii=1,nm1
          psi=psi+1.d0/ii
        enddo
        del=fact*(-log(x)+psi)
      endif
      expint=expint+del
      if(abs(del).lt.abs(expint)*EPS) return
    enddo
    stop 'series failed in expint'
  endif
  return
  END function
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
  SUBROUTINE PBEc_old(rho,t,ec_unif,drs_ec_unif,H,drs_H,dt_H)
  implicit none
  real(8),intent(IN) :: rho,t
  real(8),intent(OUT) :: ec_unif,drs_ec_unif,H,drs_H,dt_H
  real(8),parameter :: pi=3.141592653589793d0
  integer,parameter :: p=1
  real(8),parameter :: Au=0.031091d0,a1=0.21370d0,b1=7.5957d0
  real(8),parameter :: b2=3.5876d0,b3=1.6382d0,b4=0.49294d0
  real(8),parameter :: beta=0.066725d0,gamma=0.031091d0
  real(8) rs,sqrs,kf,ks,b14,exec,Ah,Ah21,Ah421,dec_Ah,dA_H

  rs=(3d0/(4d0*pi*rho))**(1d0/3d0)
  sqrs=sqrt(rs)
  kf=(3*pi**2*rho)**(1d0/3d0)
  ks=(4*kf/pi)**0.5d0
  b14=b1*sqrs+b2*rs+b3*rs*sqrs+b4*rs**(p+1)
  ec_unif=-2*Au*(1+a1*rs)*log(1+1/(2*Au*b14))
  exec=exp(-ec_unif/gamma)
  Ah=beta/gamma/(exec-1)
  Ah21=1+Ah*t**2
  Ah421=1+Ah*t**2+Ah**2*t**4
  H=gamma*log(1+beta/gamma*t**2*Ah21/Ah421)

  drs_ec_unif=-2*Au*a1*log(1+1/(2*Au*b14)) &
 &            +(1+a1*rs)*(b1/(2*sqrs)+b2+1.5d0*b3*sqrs+(p+1)*b4*rs**p)/b14**2 &
 &            /(1+1/(2*Au*b14))

  dec_Ah=beta/gamma**2*exec/(exec-1)**2
  dA_H=gamma*(beta/gamma)*t**2*(t**2/Ah421-Ah21*(t**2+2*Ah*t**4)/Ah421**2) &
 &     /(1+beta/gamma*t**2*Ah21/Ah421)
  drs_H=drs_ec_unif*dec_Ah*dA_H

  dt_H=gamma*((beta/gamma)*2*t*Ah21/Ah421 &
 &            +beta/gamma*t**2*(2*Ah*t*Ah421-Ah21*(2*Ah*t+4*Ah**2*t**3))/Ah421**2) &
 &    /(1+beta/gamma*t**2*Ah21/Ah421)
  return
  end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
! Perdew-Wang correlation energy

  SUBROUTINE PWc(rhoa,rhob,ec,dec_drhoa,dec_drhob)
  implicit none
  real(8),parameter :: Pi=3.141592653589793d0
  real(8),intent(IN) :: rhoa,rhob
  real(8),intent(OUT) :: ec,dec_drhoa,dec_drhob
  real(8),parameter :: pec0=1d0,pec1=1d0,par=1d0
  real(8),parameter :: Aec0=0.031091d0,Aec1=0.015545d0,Aar=0.016887d0
  real(8),parameter :: a1ec0=0.21370d0,a1ec1=0.20548d0,a1ar=0.11125d0
  real(8),parameter :: b1ec0=7.5957d0,b1ec1=14.1189d0,b1ar=10.357d0
  real(8),parameter :: b2ec0=3.5876d0,b2ec1=6.1977d0,b2ar=3.6231d0
  real(8),parameter :: b3ec0=1.6382d0,b3ec1=3.3662d0,b3ar=0.88026d0
  real(8),parameter :: b4ec0=0.49294d0,b4ec1=0.62517d0,b4ar=0.49671d0
  real(8),parameter :: f2d=1.709921d0,eps=1d-15
  real(8) rs,sqrs,zeta,fzeta,dfzeta,Q0ec0,Q0ec1,Q0ar,Q1ec0,Q1ec1,Q1ar
  real(8) Qdec0,Qdec1,Qdar,Gec0,Gec1,Gar,dGec0,dGec1,dGar,decdrs,decdzeta

  rs=(3d0/(4d0*Pi*(rhoa+rhob)))**(1d0/3d0)
  sqrs=sqrt(rs)
  zeta=(rhoa-rhob)/(rhoa+rhob)
  if(abs(zeta-1d0) < eps) zeta=1d0-eps
  if(abs(zeta+1d0) < eps) zeta=-1d0+eps
  fzeta=((1+zeta)**(4d0/3d0)+(1-zeta)**(4d0/3d0)-2)/(2d0**(4d0/3d0)-2d0)
  dfzeta=4d0/3d0*((1+zeta)**(1d0/3d0)-(1-zeta)**(1d0/3d0))/(2d0**(4d0/3d0)-2d0)

  Q0ec0=-2*Aec0*(1+a1ec0*rs)
  Q0ec1=-2*Aec1*(1+a1ec1*rs)
  Q0ar =-2*Aar *(1+a1ar *rs)
  Q1ec0=2*Aec0*(b1ec0*sqrs+b2ec0*rs+b3ec0*rs*sqrs+b4ec0*rs**(pec0+1))
  Q1ec1=2*Aec1*(b1ec1*sqrs+b2ec1*rs+b3ec1*rs*sqrs+b4ec1*rs**(pec1+1))
  Q1ar =2*Aar *(b1ar *sqrs+b2ar *rs+b3ar *rs*sqrs+b4ar *rs**(par +1))
  Qdec0=Aec0*(b1ec0/sqrs+2*b2ec0+3*b3ec0*sqrs+2*(pec0+1)*b4ec0*rs**pec0)
  Qdec1=Aec1*(b1ec1/sqrs+2*b2ec1+3*b3ec1*sqrs+2*(pec1+1)*b4ec1*rs**pec1)
  Qdar =Aar *(b1ar /sqrs+2*b2ar +3*b3ar *sqrs+2*(par +1)*b4ar *rs**par )
  Gec0=Q0ec0*log(1+1d0/Q1ec0)
  Gec1=Q0ec1*log(1+1d0/Q1ec1)
  Gar =Q0ar *log(1+1d0/Q1ar )
  dGec0=-2*Aec0*a1ec0*log(1+1d0/Q1ec0)-Q0ec0*Qdec0/(Q1ec0**2+Q1ec0)
  dGec1=-2*Aec1*a1ec1*log(1+1d0/Q1ec1)-Q0ec1*Qdec1/(Q1ec1**2+Q1ec1)
  dGar =-2*Aar *a1ar *log(1+1d0/Q1ar )-Q0ar *Qdar /(Q1ar **2+Q1ar )

  ec=Gec0-Gar*fzeta/f2d*(1-zeta**4)+(Gec1-Gec0)*fzeta*zeta**4
  decdrs=dGec0*(1-fzeta*zeta**4)+dGec1*fzeta*zeta**4-dGar*fzeta/f2d*(1-zeta**4)
  decdzeta=dfzeta*(-Gar/f2d*(1-zeta**4)+(Gec1-Gec0)*zeta**4) &
&         +4*zeta**3*(fzeta*(Gec1-Gec0)+Gar*fzeta/f2d)

  dec_drhoa=(-rs/3*decdrs-(zeta-1)*decdzeta)/(rhoa+rhob)
  dec_drhob=(-rs/3*decdrs-(zeta+1)*decdzeta)/(rhoa+rhob)

  return
  end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
! Perdew-Burke-Ernzerhof correlation energy
  SUBROUTINE PBEc(rhoa,rhob,grhoa,grhob,ec,dec_drhoa,dec_drhob,dec_dgrhoa,dec_dgrhob)
  implicit none
  real(8),parameter :: Pi=3.141592653589793d0
  real(8),parameter :: beta=0.066725d0,gamma=0.031091d0
  real(8),parameter :: eps=1d-12
  real(8) :: rhoa,rhob,grhoa(3),grhob(3)
  real(8) :: ec,dec_drhoa,dec_drhob,dec_dgrhoa(3),dec_dgrhob(3)
  real(8) :: ec_PW,dec_drhoa_PW,dec_drhob_PW
  real(8) :: rho,grho(3),rs,kf,ks,agrho,zeta,phi,t,exp_ec_PW,A,At2,At21,At421,bgt2A,H
  real(8) :: dH_dphi,dH_dt,dH_dA,dphi_dzeta,dzeta_drhoa,dzeta_drhob,dphi_drhoa,dphi_drhob
  real(8) :: dks_drhoa,dks_drhob,dA_drhoa,dA_drhob,dt_drhoa,dt_drhob
  real(8) :: dH_drhoa,dH_drhob,dH_dgrhoa(3),dH_dgrhob(3)

  call PWc(rhoa,rhob,ec_PW,dec_drhoa_PW,dec_drhob_PW)

  rho=rhoa+rhob
  grho(:)=grhoa(:)+grhob(:)
  rs=(3d0/(4d0*Pi*rho))**(1d0/3d0)
  kf=(3*Pi**2*rho)**(1d0/3d0)
  ks=sqrt(4*kf/Pi)
  agrho=sqrt(grho(1)**2+grho(2)**2+grho(3)**2)
  zeta=(rhoa-rhob)/rho
  if(abs(zeta-1d0) < eps) zeta=1d0-eps
  if(abs(zeta+1d0) < eps) zeta=-1d0+eps
  phi=((1+zeta)**(2d0/3d0)+(1-zeta)**(2d0/3d0))/2d0
  t=agrho/(2*phi*ks*rho)

  exp_ec_PW=exp(-ec_PW/gamma/phi**3)
  A=beta/gamma/(exp_ec_PW-1)
  At2=A*t*t
  At21=1+At2
  At421=At21+At2**2
  bgt2A=1+beta/gamma*t**2*At21/At421
  H=gamma*phi**3*log(bgt2A)

  dH_dphi=3*H/phi
  dH_dt=gamma*phi**3*beta/gamma*2*t*(1+2*At2)/At421**2/bgt2A
  dH_dA=gamma*phi**3*(-beta/gamma)*t**4*At2*(2+At2)/At421**2/bgt2A
  dphi_dzeta=(1/((1+zeta)**2+eps)**(1d0/6d0)-1/((1-zeta)**2+eps)**(1d0/6d0))/3
  dzeta_drhoa=2*rhob/rho**2
  dzeta_drhob=-2*rhoa/rho**2
  dphi_drhoa=dzeta_drhoa*dphi_dzeta
  dphi_drhob=dzeta_drhob*dphi_dzeta

  dks_drhoa=ks/(6*rho)
  dks_drhob=ks/(6*rho)
  dA_drhoa=A**2/beta*exp_ec_PW*(dec_drhoa_PW/phi**3-3*ec_PW/phi**4*dphi_drhoa)
  dA_drhob=A**2/beta*exp_ec_PW*(dec_drhob_PW/phi**3-3*ec_PW/phi**4*dphi_drhob)
  dt_drhoa=-t/phi*dphi_drhoa-t/ks*dks_drhoa-t/rho
  dt_drhob=-t/phi*dphi_drhob-t/ks*dks_drhob-t/rho
  dH_drhoa=dphi_drhoa*dH_dphi+dt_drhoa*dH_dt+dA_drhoa*dH_dA
  dH_drhob=dphi_drhob*dH_dphi+dt_drhob*dH_dt+dA_drhob*dH_dA

  dH_dgrhoa(:)=dH_dt*t*grho(:)/agrho**2
  dH_dgrhob(:)=dH_dgrhoa(:)

  ec=ec_PW+H
  dec_drhoa=dec_drhoa_PW+dH_drhoa
  dec_drhob=dec_drhob_PW+dH_drhob
  dec_dgrhoa(:)=dH_dgrhoa(:)
  dec_dgrhob(:)=dH_dgrhob(:)

  return
  end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
  use Global_Variables
  use communication
  implicit none
  integer,intent(in)    :: GS_RT
  real(8),intent(inout) :: rho_s(NL),tau_s(NL),j_s(NL,3),grho_s(NL,3),lrho_s(NL)
  integer :: ikb,ik,ib,i
  real(8) :: tau_s_l(NL),j_s_l(NL,3),ss(3)
  complex(8) :: zs(3)
  integer :: thr_id,omp_get_thread_num

!  allocate(tau_s_l_omp(NL,0:NUMBER_THREADS-1),j_s_l_omp(NL,3,0:NUMBER_THREADS-1)) ! sato
  rho_s=rho*0.5d0
  if(flag_nlcc)rho_s = rho_s + 0.5d0*rho_nlcc

  tau_s_l_omp=0d0
  j_s_l_omp=0d0

  if(GS_RT == calc_mode_gs)then

    select case(Nd)
    case(4)

      thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ik,ib,zs,i)
      do ikb=1,NKB
        ik=ik_table(ikb) ; ib=ib_table(ikb)  
        do i=1,NL
          zs(1)=nabx(1)*(zu_GS(ifdx(1,i),ib,ik)-zu_GS(ifdx(-1,i),ib,ik))&
            &  +nabx(2)*(zu_GS(ifdx(2,i),ib,ik)-zu_GS(ifdx(-2,i),ib,ik))&
            &  +nabx(3)*(zu_GS(ifdx(3,i),ib,ik)-zu_GS(ifdx(-3,i),ib,ik))&
            &  +nabx(4)*(zu_GS(ifdx(4,i),ib,ik)-zu_GS(ifdx(-4,i),ib,ik))&
            &  +zI*kAc0(ik,1)*zu_GS(i,ib,ik)
          zs(2)=naby(1)*(zu_GS(ifdy(1,i),ib,ik)-zu_GS(ifdy(-1,i),ib,ik))&
            &  +naby(2)*(zu_GS(ifdy(2,i),ib,ik)-zu_GS(ifdy(-2,i),ib,ik))&
            &  +naby(3)*(zu_GS(ifdy(3,i),ib,ik)-zu_GS(ifdy(-3,i),ib,ik))&
            &  +naby(4)*(zu_GS(ifdy(4,i),ib,ik)-zu_GS(ifdy(-4,i),ib,ik))&
            &  +zI*kAc0(ik,2)*zu_GS(i,ib,ik)            
          zs(3)=nabz(1)*(zu_GS(ifdz(1,i),ib,ik)-zu_GS(ifdz(-1,i),ib,ik))&
            &  +nabz(2)*(zu_GS(ifdz(2,i),ib,ik)-zu_GS(ifdz(-2,i),ib,ik))&
            &  +nabz(3)*(zu_GS(ifdz(3,i),ib,ik)-zu_GS(ifdz(-3,i),ib,ik))&
            &  +nabz(4)*(zu_GS(ifdz(4,i),ib,ik)-zu_GS(ifdz(-4,i),ib,ik))&
            &  +zI*kAc0(ik,3)*zu_GS(i,ib,ik)   
          tau_s_l_omp(i,thr_id)=tau_s_l_omp(i,thr_id) &
            &+(abs(zs(1))**2+abs(zs(2))**2+abs(zs(3))**2)*(occ(ib,ik)*0.5d0)*0.5d0
          j_s_l_omp(i,1,thr_id)=j_s_l_omp(i,1,thr_id)+imag(conjg(zu_GS(i,ib,ik))*zs(1))*(occ(ib,ik)*0.5d0)
          j_s_l_omp(i,2,thr_id)=j_s_l_omp(i,2,thr_id)+imag(conjg(zu_GS(i,ib,ik))*zs(2))*(occ(ib,ik)*0.5d0)
          j_s_l_omp(i,3,thr_id)=j_s_l_omp(i,3,thr_id)+imag(conjg(zu_GS(i,ib,ik))*zs(3))*(occ(ib,ik)*0.5d0)
        enddo
      end do
!$omp end parallel

    case default
      call err_finalize('Nd /= 4')
    end select
  
  else  if(GS_RT == calc_mode_rt)then

    select case(Nd)
    case(4)

      thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ik,ib,zs,i)
      do ikb=1,NKB
        ik=ik_table(ikb) ; ib=ib_table(ikb)  
        do i=1,NL
          zs(1)=nabx(1)*(zu(ifdx(1,i),ib,ik)-zu(ifdx(-1,i),ib,ik))&
            &  +nabx(2)*(zu(ifdx(2,i),ib,ik)-zu(ifdx(-2,i),ib,ik))&
            &  +nabx(3)*(zu(ifdx(3,i),ib,ik)-zu(ifdx(-3,i),ib,ik))&
            &  +nabx(4)*(zu(ifdx(4,i),ib,ik)-zu(ifdx(-4,i),ib,ik))&
            &  +zI*kAc0(ik,1)*zu(i,ib,ik)
          zs(2)=naby(1)*(zu(ifdy(1,i),ib,ik)-zu(ifdy(-1,i),ib,ik))&
            &  +naby(2)*(zu(ifdy(2,i),ib,ik)-zu(ifdy(-2,i),ib,ik))&
            &  +naby(3)*(zu(ifdy(3,i),ib,ik)-zu(ifdy(-3,i),ib,ik))&
            &  +naby(4)*(zu(ifdy(4,i),ib,ik)-zu(ifdy(-4,i),ib,ik))&
            &  +zI*kAc0(ik,2)*zu(i,ib,ik)
          zs(3)=nabz(1)*(zu(ifdz(1,i),ib,ik)-zu(ifdz(-1,i),ib,ik))&
            &  +nabz(2)*(zu(ifdz(2,i),ib,ik)-zu(ifdz(-2,i),ib,ik))&
            &  +nabz(3)*(zu(ifdz(3,i),ib,ik)-zu(ifdz(-3,i),ib,ik))&
            &  +nabz(4)*(zu(ifdz(4,i),ib,ik)-zu(ifdz(-4,i),ib,ik))&
            &  +zI*kAc0(ik,3)*zu(i,ib,ik)
          tau_s_l_omp(i,thr_id)=tau_s_l_omp(i,thr_id) &
            &+(abs(zs(1))**2+abs(zs(2))**2+abs(zs(3))**2)*(occ(ib,ik)*0.5d0)*0.5d0
          j_s_l_omp(i,1,thr_id)=j_s_l_omp(i,1,thr_id)+imag(conjg(zu(i,ib,ik))*zs(1))*(occ(ib,ik)*0.5d0)
          j_s_l_omp(i,2,thr_id)=j_s_l_omp(i,2,thr_id)+imag(conjg(zu(i,ib,ik))*zs(2))*(occ(ib,ik)*0.5d0)
          j_s_l_omp(i,3,thr_id)=j_s_l_omp(i,3,thr_id)+imag(conjg(zu(i,ib,ik))*zs(3))*(occ(ib,ik)*0.5d0)
        enddo
      end do
!$omp end parallel

    case default
      call err_finalize('Nd /= 4')
    end select

  else
    call err_finalize('error in meta GGA')
  end if

  tau_s_l(:) = tau_s_l_omp(:,0)
  j_s_l(:,1) = j_s_l_omp(:,1,0)
  j_s_l(:,2) = j_s_l_omp(:,2,0)
  j_s_l(:,3) = j_s_l_omp(:,3,0)

  do thr_id = 1, NUMBER_THREADS-1
!$omp parallel do    
    do i=1,NL
      tau_s_l(i) = tau_s_l(i) + tau_s_l_omp(i,thr_id)
      j_s_l(i,1) = j_s_l(i,1) + j_s_l_omp(i,1,thr_id)
      j_s_l(i,2) = j_s_l(i,2) + j_s_l_omp(i,2,thr_id)
      j_s_l(i,3) = j_s_l(i,3) + j_s_l_omp(i,3,thr_id)
    end do
  end do

  call comm_summation(tau_s_l,tau_s,NL,proc_group(2))
  call comm_summation(j_s_l,j_s,NL*3,proc_group(2))

  if(flag_nlcc)tau_s = tau_s + 0.5d0*tau_nlcc

  select case(Nd)
  case(4)
!$omp parallel do  private(ss)  
    do i=1,NL 
      grho_s(i,1)=nabx(1)*(rho_s(ifdx(1,i))-rho_s(ifdx(-1,i))) &
        &        +nabx(2)*(rho_s(ifdx(2,i))-rho_s(ifdx(-2,i))) &
        &        +nabx(3)*(rho_s(ifdx(3,i))-rho_s(ifdx(-3,i))) &
        &        +nabx(4)*(rho_s(ifdx(4,i))-rho_s(ifdx(-4,i)))
      grho_s(i,2)=naby(1)*(rho_s(ifdy(1,i))-rho_s(ifdy(-1,i))) &
        &        +naby(2)*(rho_s(ifdy(2,i))-rho_s(ifdy(-2,i))) &
        &        +naby(3)*(rho_s(ifdy(3,i))-rho_s(ifdy(-3,i))) &
        &        +naby(4)*(rho_s(ifdy(4,i))-rho_s(ifdy(-4,i)))
      grho_s(i,3)=nabz(1)*(rho_s(ifdz(1,i))-rho_s(ifdz(-1,i))) &
        &        +nabz(2)*(rho_s(ifdz(2,i))-rho_s(ifdz(-2,i))) &
        &        +nabz(3)*(rho_s(ifdz(3,i))-rho_s(ifdz(-3,i))) &
        &        +nabz(4)*(rho_s(ifdz(4,i))-rho_s(ifdz(-4,i)))
      ss(1)=lapx(0)*rho_s(i)&
          &+lapx(1)*(rho_s(ifdx(1,i))+rho_s(ifdx(-1,i))) &
          &+lapx(2)*(rho_s(ifdx(2,i))+rho_s(ifdx(-2,i))) &
          &+lapx(3)*(rho_s(ifdx(3,i))+rho_s(ifdx(-3,i))) &
          &+lapx(4)*(rho_s(ifdx(4,i))+rho_s(ifdx(-4,i)))
      ss(2)=lapy(0)*rho_s(i)&
          &+lapy(1)*(rho_s(ifdy(1,i))+rho_s(ifdy(-1,i))) &
          &+lapy(2)*(rho_s(ifdy(2,i))+rho_s(ifdy(-2,i))) &
          &+lapy(3)*(rho_s(ifdy(3,i))+rho_s(ifdy(-3,i))) &
          &+lapy(4)*(rho_s(ifdy(4,i))+rho_s(ifdy(-4,i)))
      ss(3)=lapz(0)*rho_s(i)&
          &+lapz(1)*(rho_s(ifdz(1,i))+rho_s(ifdz(-1,i))) &
          &+lapz(2)*(rho_s(ifdz(2,i))+rho_s(ifdz(-2,i))) &
          &+lapz(3)*(rho_s(ifdz(3,i))+rho_s(ifdz(-3,i))) &
          &+lapz(4)*(rho_s(ifdz(4,i))+rho_s(ifdz(-4,i)))
      lrho_s(i)=ss(1)+ss(2)+ss(3)
    enddo

  case default
    call err_finalize('Nd /= 4')
  end select

!sato
! Symmetry
  select case(crystal_structure)
  case("diamond")
     if(Sym == 4)then
        tau_s(:)=tau_s(:)*0.25d0
        j_s(:,1:3)=j_s(:,1:3)*0.25d0

!tau_s_l(NL),j_s_l(NL,3),ss(3)
! 1.T_3
!$omp parallel do  private(ss,i)  
        do i=1,NL
           tau_s_l(i)=tau_s(i)+tau_s(itable_sym(1,i))
           j_s_l(i,1)=j_s(i,1)-j_s(itable_sym(1,i),2)
           j_s_l(i,2)=j_s(i,2)+j_s(itable_sym(1,i),1)
           j_s_l(i,3)=j_s(i,3)+j_s(itable_sym(1,i),3)
        end do
! 2.T3*T3
!$omp parallel do  private(ss,i)  
        do i=1,NL
           tau_s(i)=tau_s_l(i)+tau_s_l(itable_sym(2,i))
           j_s(i,1)=j_s_l(i,1)-j_s_l(itable_sym(2,i),1)
           j_s(i,2)=j_s_l(i,2)-j_s_l(itable_sym(2,i),2)
           j_s(i,3)=j_s_l(i,3)+j_s_l(itable_sym(2,i),3)
        end do

     else if(Sym == 8)then
        tau_s(:)=tau_s(:)/32d0
        j_s(:,1:3)=j_s(:,1:3)/32d0

! 1.T_4
!$omp parallel do  private(i)  
        do i=1,NL
           tau_s_l(i)=tau_s(i)+tau_s(itable_sym(4,i))
           j_s_l(i,1)=j_s(i,1)+j_s(itable_sym(4,i),1)
           j_s_l(i,2)=j_s(i,2)+j_s(itable_sym(4,i),2)
           j_s_l(i,3)=j_s(i,3)+j_s(itable_sym(4,i),3)
        end do        

! 2.T_3*T_3
!$omp parallel do  private(i)  
        do i=1,NL
           tau_s(i)=tau_s_l(i)+tau_s_l(itable_sym(5,i))
           j_s(i,1)=j_s_l(i,1)-j_s_l(itable_sym(5,i),1)
           j_s(i,2)=j_s_l(i,2)-j_s_l(itable_sym(5,i),2)
           j_s(i,3)=j_s_l(i,3)+j_s_l(itable_sym(5,i),3)
        end do        
        
! 3.T_3
!$omp parallel do  private(i)  
        do i=1,NL
           tau_s_l(i)=tau_s(i)+tau_s(itable_sym(3,i))
           j_s_l(i,1)=j_s(i,1)-j_s(itable_sym(3,i),2)
           j_s_l(i,2)=j_s(i,2)+j_s(itable_sym(3,i),1)
           j_s_l(i,3)=j_s(i,3)+j_s(itable_sym(3,i),3)
        end do

! 4.T_1
!$omp parallel do  private(i)  
        do i=1,NL
           tau_s(i)=tau_s_l(i)+tau_s_l(itable_sym(1,i))
           j_s(i,1)=j_s_l(i,1)+j_s_l(itable_sym(1,i),2)
           j_s(i,2)=j_s_l(i,2)+j_s_l(itable_sym(1,i),1)
           j_s(i,3)=j_s_l(i,3)+j_s_l(itable_sym(1,i),3)
        end do        

! 5.T_2
!$omp parallel do  private(i)  
        do i=1,NL
           tau_s_l(i)=tau_s(i)+tau_s(itable_sym(2,i))
           j_s_l(i,1)=j_s(i,1)-j_s(itable_sym(2,i),2)
           j_s_l(i,2)=j_s(i,2)-j_s(itable_sym(2,i),1)
           j_s_l(i,3)=j_s(i,3)+j_s(itable_sym(2,i),3)
        end do        

! 6. I
!$omp parallel do  private(i)  
        do i=1,NL
           tau_s(i)=tau_s_l(i)
           j_s(i,1:3)=j_s_l(i,1:3)
        end do        
        
     else if(Sym /= 1)then
        call err_finalize('Bad crystal structure')
     end if

  case("diamond2")
     if(Sym == 8)then
        tau_s(:)=tau_s(:)*6.25d-2
        j_s(:,1:3)=j_s(:,1:3)*6.25d-2

!tau_s_l(NL),j_s_l(NL,3),ss(3)
!$omp parallel do  private(i)  
        do i=1,NL
           tau_s_l(i)=tau_s(i)+tau_s(itable_sym(4,i))
           j_s_l(i,1)=j_s(i,1)+j_s(itable_sym(4,i),2)
           j_s_l(i,2)=j_s(i,2)-j_s(itable_sym(4,i),1)
           j_s_l(i,3)=j_s(i,3)+j_s(itable_sym(4,i),3)
        end do        

!$omp parallel do  private(i)  
        do i=1,NL
           tau_s(i)=tau_s_l(i)+tau_s_l(itable_sym(2,i))
           j_s(i,1)=j_s_l(i,1)+j_s_l(itable_sym(2,i),1)
           j_s(i,2)=j_s_l(i,2)-j_s_l(itable_sym(2,i),2)
           j_s(i,3)=j_s_l(i,3)+j_s_l(itable_sym(2,i),3)
        end do        
        
!$omp parallel do  private(i)  
        do i=1,NL
           tau_s_l(i)=tau_s(i)+tau_s(itable_sym(3,i))
           j_s_l(i,1)=j_s(i,1)+j_s(itable_sym(3,i),1)
           j_s_l(i,2)=j_s(i,2)+j_s(itable_sym(3,i),2)
           j_s_l(i,3)=j_s(i,3)+j_s(itable_sym(3,i),3)
        end do

!$omp parallel do  private(i)  
        do i=1,NL
           tau_s(i)=tau_s_l(i)+tau_s_l(itable_sym(1,i))
           j_s(i,1)=j_s_l(i,1)-j_s_l(itable_sym(1,i),1)
           j_s(i,2)=j_s_l(i,2)+j_s_l(itable_sym(1,i),2)
           j_s(i,3)=j_s_l(i,3)+j_s_l(itable_sym(1,i),3)
        end do        

     else if(Sym /= 1)then
        call err_finalize('Bad crystal structure')
     end if
  case default
     if(Sym /= 1)call err_finalize('Bad crystal structure')
  end select
!sato

  return
End Subroutine rho_j_tau

end subroutine Exc_Cor
