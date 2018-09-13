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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine calc_vloc(pp,dvloc_g,gx,gy,gz,ng,ng_s,ng_e,ngzero)
  use salmon_global,only : nelem
  use salmon_pp,only : pp_info
  implicit none
  real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
  type(pp_info) :: pp
  integer,intent(in) :: ng,ng_s,ng_e
  integer,intent(in) :: ngzero
  real(8),intent(in) :: gx(ng),gy(ng),gz(ng)
  complex(8),intent(out) :: dvloc_g(ng_s:ng_e,nelem)
  integer :: i,ik,n
  real(8) :: dr
  real(8) :: g2sq
  real(8) :: vloc_av
  real(8) :: s,r

!$omp parallel
!$omp do private(ik,n,g2sq,s,r,dr,i,vloc_av) collapse(2)
  do ik=1,nelem
    do n=ng_s,ng_e
      g2sq=sqrt(gx(n)**2+gy(n)**2+gz(n)**2)
      s=0.d0
      if (n == ngzero) then
        do i=2,pp%nrloc(ik)
          r=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
          dr=pp%rad(i,ik)-pp%rad(i-1,ik)
          vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
          s=s+4d0*pi*(r**2*vloc_av+r*pp%zps(ik))*dr
        enddo
      else
        do i=2,pp%nrloc(ik)
          r=0.5d0*(pp%rad(i,ik)+pp%rad(i-1,ik))
          dr=pp%rad(i,ik)-pp%rad(i-1,ik)
          vloc_av = 0.5d0*(pp%vloctbl(i,ik)+pp%vloctbl(i-1,ik))
          s=s+4d0*pi*sin(g2sq*r)/g2sq*(r*vloc_av+pp%zps(ik))*dr !Vloc - coulomb
        enddo
      endif
      dvloc_g(n,ik)=s
    enddo
  enddo
!$omp end do
!$omp end parallel

end subroutine calc_vloc

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130


subroutine calc_vpsl(pp,rhoion_g,vpsl_ia,vpsl,dvloc_g,  &
                     ngzero,gx,gy,gz,ng,ng_s,ng_e,nl,alxyz,lx,ly,lz,hx,hy,hz)
  use salmon_global,only : natom, nelem, kion, rion
  use salmon_parallel,only : nproc_group_tdks
  use salmon_communication, only: comm_summation
  use salmon_pp,only : pp_info
  implicit none
  real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
  type(pp_info) :: pp
  integer,intent(in) :: ngzero
  integer,intent(in) :: ng,ng_s,ng_e
  complex(8),intent(in) :: dvloc_g(ng_s:ng_e,nelem)
  real(8),intent(in) :: gx(ng),gy(ng),gz(ng)
  integer,intent(in) :: nl
  real(8),intent(in) :: alxyz
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  real(8),intent(in) :: hx,hy,hz
  complex(8),intent(out) :: rhoion_g(ng_s:ng_e)
  real(8),intent(out) :: vpsl_ia(nl,natom)
  real(8),intent(out) :: vpsl(nl)
  integer :: a,i,n,ik
  real(8) :: gd
  real(8) :: g2
  complex(8) :: vion_g_ia(ng_s:ng_e,natom),tmp_exp !, Vion_G(NG_s:NG_e)
  real(8) :: vpsl_ia_l(nl,natom)
  real(8) :: gr
  complex(8),parameter :: zi=(0.d0,1.d0)

  !(Local pseudopotential: Vlocal in G-space(=Vion_G))
  vion_g_ia=0.d0
 !Vion_G   =0.d0
  rhoion_g =0.d0
!$omp parallel private(a,ik)
  do a=1,natom
    ik=kion(a)
!$omp do private(n,g2,gd,tmp_exp)
    do n=ng_s,ng_e
      gd=gx(n)*rion(1,a)+gy(n)*rion(2,a)+gz(n)*rion(3,a)
      tmp_exp = exp(-zi*gd)/alxyz
     !Vion_G(n)     = Vion_G(n)      + dvloc_g(n,ik)*tmp_exp
      vion_g_ia(n,a)= vion_g_ia(n,a) + dvloc_g(n,ik)*tmp_exp
      rhoion_g(n)   = rhoion_g(n) + pp%zps(ik)*tmp_exp
      if(n == ngzero) cycle
      !(add coulomb as dvloc_g is given by Vloc - coulomb)
      g2=gx(n)**2+gy(n)**2+gz(n)**2
     !Vion_G(n)     = Vion_G(n)      -4d0*pi/g2*pp%zps(ik)*tmp_exp
      vion_g_ia(n,a)= vion_g_ia(n,a) -4d0*pi/g2*pp%zps(ik)*tmp_exp
    enddo
!$omp end do
  enddo
!$omp end parallel

  !(Local pseudopotential: Vlocal(=Vpsl) in real-space)
  vpsl_ia_l=0.d0
 !Vpsl_l   =0.d0
!$omp parallel private(n)
  do n=ng_s,ng_e
!$omp do private(i,gr,a,tmp_exp)
    do i=1,NL
      gr = gx(n)*lx(i)*hx+gy(n)*ly(i)*hy+gz(n)*lz(i)*hz
      tmp_exp = exp(zi*gr)
      !vpsl_l(i) = vpsl_l(i) + vion_G(n)*tmp_exp
      do a=1,natom
        vpsl_ia_l(i,a)= vpsl_ia_l(i,a)+ vion_g_ia(n,a)*tmp_exp
      enddo
    enddo
!$omp end do
  enddo
!$omp end parallel

 !call comm_summation(vpsl_l,vpsl,nl,nproc_group_tdks)
  call comm_summation(vpsl_ia_l,vpsl_ia,nl*natom,nproc_group_tdks)

!$omp parallel
!$omp do private(i)
  do i=1,nl
    vpsl(i) = sum(vpsl_ia(i,:))
  enddo
!$omp end do
!$omp end parallel



end subroutine calc_vpsl
