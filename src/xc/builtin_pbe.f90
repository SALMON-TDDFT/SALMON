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
module builtin_pbe
  use salmon_math, only: erfc_salmon
  implicit none
  
  real(8),parameter :: Pi=3.141592653589793d0
  
contains
  
  subroutine exc_cor_pbe(nl, rho, grho_s, exc, eexc, vexc, nd, ifdx, ifdy, ifdz, nabx, naby, nabz)
    implicit none
    integer, intent(in) :: nl
    real(8), intent(in) :: rho(nl), grho_s(nl, 3)
    real(8), intent(out) :: exc(nl), eexc(nl), vexc(nl)
    integer, intent(in) :: nd
    integer, intent(in) :: ifdx(-nd:nd, 1:nl)
    integer, intent(in) :: ifdy(-nd:nd, 1:nl)
    integer, intent(in) :: ifdz(-nd:nd, 1:nl)
    real(8), intent(in) :: nabx(-nd:nd)
    real(8), intent(in) :: naby(-nd:nd)
    real(8), intent(in) :: nabz(-nd:nd)
    
    real(8) :: agrho_s(NL)
    integer :: i
    real(8) :: omega
    real(8) :: FxHSE,drho_FxHSE,ds_FxHSE,ec_unif,drs_ec_unif,H,drs_H,dt_H
    real(8) :: rho_t,rs,kf,ks,grho(NL,3),agrho(NL),srho,trho
    real(8) :: exlda,Ex(NL),Ec(NL),Vx(NL),Vc(NL)
    real(8) :: dFxdrho_tbl(NL,3),dHksrho_tbl(NL,3)
    

    ! call rho_j_tau(GS_RT,rho_s,tau_s,j_s,grho_s,lrho_s)
    
    ! rho -> rho_s
    ! |tau_s_l_omp, j_s_l_omp, thr_id,
    ! ?ND branch and loop:
    ! | (nabla + ikAc) zu_GS -> zs
    ! | zs**2 * occ -> tau_s_l_omp
    ! | zu_GS.T * 1j * zs * occ -> j_s_l_omp
    ! tau_s_l_omp -> tau_s_l [reduction omp]
    ! j_s_l_omp -> j_s_l [reduction omp]
    ! tau_s_l -> tau_s [reduction mpi]
    ! j_s_l -> j_s [reduction mpi]
    ! tau_s += 0.5d0*tau_nlcc ???
    ! ?Sym branch and loop:
    ! | nabla * rho_s -> grho_s
    ! | lap * rho_s -> ss
    ! | ss**2 -> lrho_s
    ! Symmetry Expansion!
    
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
    Exc=(Ex+Ec)
    Eexc=rho*(Ex+Ec)
    Vexc=Vx+Vc
    return
  End Subroutine Exc_Cor_PBE


  SUBROUTINE HSEFx(rho,ss,omega,FxHSE,drho_FxHSE,ds_FxHSE)
  implicit none

  real(8),intent(IN) :: rho,ss,omega
  real(8),intent(OUT) :: FxHSE,drho_FxHSE,ds_FxHSE

  !real(8),parameter :: pi=3.1415926535897932d0
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
 &     *erfc_salmon(1.5d0*s*sqrt(H/A))
    G_b=(15d0/16d0*sqpi*s2)/DHs72
    EG=-(0.75d0*pi+G_a)/G_b
    ds_EGs=-14d0/5d0*sqpi*DHs52*ds_Hs &
 &         +sqpi*sqrt(A)/5*DHs52  &
 &         *exp(9*H*s2/(4*A))*erfc_salmon(1.5d0*s*sqrt(H/A)) &
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
      exer=pi*exp(HA94)*erfc_salmon(sqrt(HA94))
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
 &   -3d0/4d0*sqrt(A)*pi*exp(A94*(H*s2+w2))*erfc_salmon(sqrt(A94*(H*s2+w2))) &
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
  
  
  SUBROUTINE PBEc_old(rho,t,ec_unif,drs_ec_unif,H,drs_H,dt_H)
  implicit none
  real(8),intent(IN) :: rho,t
  real(8),intent(OUT) :: ec_unif,drs_ec_unif,H,drs_H,dt_H
  !real(8),parameter :: pi=3.141592653589793d0
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


  
end module builtin_pbe
