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

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_mps(ppg)
  use salmon_global,only : natom
  use salmon_pp,only : pp_grid
  implicit none 
  type(pp_grid) :: ppg

  allocate(ppg%mps(natom))

end subroutine init_mps
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_jxyz(ppg)
  use salmon_global,only : natom
  use salmon_pp,only : pp_grid
  implicit none 
  type(pp_grid) :: ppg

  allocate(ppg%jxyz(3,ppg%nps,natom))
  allocate(ppg%jxx( ppg%nps,natom))
  allocate(ppg%jyy( ppg%nps,natom))
  allocate(ppg%jzz( ppg%nps,natom))

end subroutine init_jxyz
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine finalize_jxyz(ppg)
  use salmon_pp,only : pp_grid
  implicit none 
  type(pp_grid) :: ppg

  deallocate(ppg%jxyz)
  deallocate(ppg%jxx,ppg%jyy,ppg%jzz)

end subroutine finalize_jxyz
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

subroutine calc_mps(pp,ppg,alx,aly,alz,lx,ly,lz,nl,mx,my,mz,ml,hx,hy,hz)
  use salmon_global,only : natom,kion,rion,iperiodic,domain_parallel
  use salmon_pp,only : pp_info,pp_grid
  implicit none
  type(pp_info) :: pp
  type(pp_grid) :: ppg
  real(8),intent(in) :: alx,aly,alz
  integer,intent(in) :: nl,ml
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  integer,intent(in) :: mx(ml),my(ml),mz(ml)
  real(8),intent(in) :: hx,hy,hz
  integer :: a,i,ik,ix,iy,iz,j
  integer :: nc
  real(8) :: tmpx,tmpy,tmpz
  real(8) :: x,y,z,r
  real(8) :: rshift(3)

  if(iperiodic==0)then
    nc=0
  else if(iperiodic==3)then
    nc=2
  end if

  if(iperiodic==0)then
    if(mod(lx(nl)-lx(1)+1,2)==1)then
      rshift(1)=0.d0
    else
      rshift(1)=-0.5d0*Hx
    end if
    if(mod(ly(nl)-ly(1)+1,2)==1)then
      rshift(2)=0.d0
    else
      rshift(2)=-0.5d0*Hy
    end if
    if(mod(lz(nl)-lz(1)+1,2)==1)then
      rshift(3)=0.d0
    else
      rshift(3)=-0.5d0*Hz
    end if
  else if(iperiodic==3)then
    if(domain_parallel=='y')then
      rshift(1)=-Hx
      rshift(2)=-Hy
      rshift(3)=-Hz
    else
      rshift(:)=0.d0
    end if
  end if


!$omp parallel
!$omp do private(a,ik,j,i,ix,iy,iz,x,y,z,r,tmpx,tmpy,tmpz)
  do a=1,natom
    ik=kion(a)
    j=0
    do ix=-nc,nc
    do iy=-nc,nc
    do iz=-nc,nc
      tmpx = rion(1,a)+ix*alx
      tmpy = rion(2,a)+iy*aly
      tmpz = rion(3,a)+iz*alz
      do i=1,ml
        x=mx(i)*Hx+rshift(1)-tmpx
        y=my(i)*Hy+rshift(2)-tmpy
        z=mz(i)*Hz+rshift(3)-tmpz
        r=sqrt(x*x+y*y+z*z)
        if (r<pp%rps(ik)) j=j+1
      enddo
    enddo
    enddo
    enddo
    ppg%mps(a)=j
  end do
!$omp end do
!$omp end parallel

  ppg%nps=maxval(ppg%mps(:))

end subroutine calc_mps

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine calc_jxyz(pp,ppg,alx,aly,alz,lx,ly,lz,nl,mx,my,mz,ml,hx,hy,hz)
  use salmon_global,only : natom,kion,rion,iperiodic,domain_parallel
  use salmon_pp,only : pp_info,pp_grid
  implicit none
  type(pp_info) :: pp
  type(pp_grid) :: ppg
  real(8),intent(in) :: alx,aly,alz
  integer,intent(in) :: nl,ml
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  integer,intent(in) :: mx(ml),my(ml),mz(ml)
  real(8),intent(in) :: hx,hy,hz
  integer :: a,i,ik,ix,iy,iz,j
  integer :: nc
  real(8) :: tmpx,tmpy,tmpz
  real(8) :: r,x,y,z
  real(8) :: rshift(3)

  if(iperiodic==0)then
    nc=0
  else if(iperiodic==3)then
    nc=2
  end if

  if(iperiodic==0)then
    if(mod(lx(nl)-lx(1)+1,2)==1)then
      rshift(1)=0.d0
    else
      rshift(1)=-0.5d0*Hx
    end if
    if(mod(ly(nl)-ly(1)+1,2)==1)then
      rshift(2)=0.d0
    else
      rshift(2)=-0.5d0*Hy
    end if
    if(mod(lz(nl)-lz(1)+1,2)==1)then
      rshift(3)=0.d0
    else
      rshift(3)=-0.5d0*Hz
    end if
  else if(iperiodic==3)then 
    if(domain_parallel=='y')then
      rshift(1)=-Hx
      rshift(2)=-Hy
      rshift(3)=-Hz
    else
      rshift(:)=0.d0
    end if
  end if

!$omp parallel
!$omp do private(a,ik,j,ix,iy,iz,tmpx,tmpy,tmpz,i,x,y,z,r)
  do a=1,natom
    ik=kion(a)
    j=0
    do ix=-nc,nc
    do iy=-nc,nc
    do iz=-nc,nc
      tmpx = rion(1,a)+ix*aLx
      tmpy = rion(2,a)+iy*aLy
      tmpz = rion(3,a)+iz*aLz
      do i=1,ml
        x=mx(i)*Hx+rshift(1)-tmpx
        y=my(i)*Hy+rshift(2)-tmpy
        z=mz(i)*Hz+rshift(3)-tmpz
        r=sqrt(x*x+y*y+z*z)
        if (r<pp%rps(ik)) then
          j=j+1
          if (j<=ppg%nps) then
            ppg%jxyz(1,j,a)=mx(i)
            ppg%jxyz(2,j,a)=my(i)
            ppg%jxyz(3,j,a)=mz(i)
            ppg%jxx( j,a)=ix
            ppg%jyy( j,a)=iy
            ppg%jzz( j,a)=iz
          endif
        endif
      enddo
    enddo
    enddo
    enddo
  end do
!$omp end do
!$omp end parallel

end subroutine calc_jxyz
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_lma_tbl(pp,ppg)
  use salmon_global,only : natom
  use salmon_pp,only : pp_info,pp_grid
  implicit none 
  type(pp_info) :: pp
  type(pp_grid) :: ppg

  allocate(ppg%lma_tbl((pp%lmax+1)**2,natom))

end subroutine init_lma_tbl
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine finalize_lma_tbl(ppg)
  use salmon_pp,only : pp_grid
  implicit none 
  type(pp_grid) :: ppg

  deallocate(ppg%lma_tbl)

end subroutine finalize_lma_tbl
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_uv(pp,ppg)
  use salmon_global,only : natom
  use salmon_pp,only : pp_info,pp_grid
  implicit none 
  type(pp_info) :: pp
  type(pp_grid) :: ppg

  allocate(ppg%ia_tbl((pp%lmax+1)**2*natom))
  allocate(ppg%rinv_uvu((pp%lmax+1)**2*natom))
  allocate(ppg%uv(ppg%nps,ppg%nlma),ppg%duv(ppg%nps,ppg%nlma,3))

end subroutine init_uv
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine finalize_uv(ppg)
  use salmon_pp,only : pp_grid
  implicit none 
  type(pp_grid) :: ppg

  deallocate(ppg%ia_tbl,ppg%rinv_uvu)
  deallocate(ppg%uv,ppg%duv)

end subroutine finalize_uv
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine set_nlma(pp,ppg)
  use salmon_global,only : natom,kion
  use salmon_pp,only : pp_info,pp_grid
  implicit none 
  type(pp_info) :: pp
  type(pp_grid) :: ppg
  integer :: lma
  integer :: a,ik,m,l

  lma=0
  do a=1,natom
    ik=kion(a)
    do l=0,pp%mlps(ik)
      if(pp%inorm(l,ik)==0) cycle
      do m=-l,l
        lma=lma+1
      enddo
    enddo
  enddo
  ppg%nlma=lma

end subroutine set_nlma
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine set_lma_tbl(pp,ppg)
  use salmon_global,only : natom,kion
  use salmon_pp,only : pp_info,pp_grid
  implicit none 
  type(pp_info) :: pp
  type(pp_grid) :: ppg
  integer :: lm,lma
  integer :: a,ik,m,l

  lma=0
  do a=1,natom
    ik=kion(a)
    lm=0
    do l=0,pp%mlps(ik)
      if(pp%inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        lma=lma+1
        ppg%lma_tbl(lm,a)=lma
        ppg%ia_tbl(lma)=a
      enddo
    enddo
  enddo

end subroutine set_lma_tbl
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine calc_uv(pp,ppg,save_udvtbl_a,save_udvtbl_b,save_udvtbl_c,save_udvtbl_d, &
                   lx,ly,lz,nl,hx,hy,hz,alx,aly,alz,  &
                   flag_use_grad_wf_on_force,property)
  use salmon_global,only : natom,kion,rion,iperiodic,domain_parallel
  use salmon_pp,only : pp_info,pp_grid
  implicit none
  real(8),parameter :: pi=3.141592653589793d0 ! copied from salmon_math
  type(pp_info) :: pp
  type(pp_grid) :: ppg
  integer,intent(in) :: nl
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  real(8),intent(in) :: hx,hy,hz
  real(8),intent(in) :: alx,aly,alz
  logical,intent(in) :: flag_use_grad_wf_on_force
  character(17),intent(in) :: property
  real(8),intent(out) :: save_udvtbl_a(pp%nrmax,0:pp%lmax,natom)
  real(8),intent(out) :: save_udvtbl_b(pp%nrmax,0:pp%lmax,natom)
  real(8),intent(out) :: save_udvtbl_c(pp%nrmax,0:pp%lmax,natom)
  real(8),intent(out) :: save_udvtbl_d(pp%nrmax,0:pp%lmax,natom)
  integer :: a,ik,j,l,lm,m
  integer :: ilma,intr,ir,lma
  real(8),allocatable :: xn(:),yn(:),an(:),bn(:),cn(:),dn(:)  
  real(8) :: dudvtbl_a(pp%nrmax,0:pp%lmax), dudvtbl_b(pp%nrmax,0:pp%lmax)
  real(8) :: dudvtbl_c(pp%nrmax,0:pp%lmax), dudvtbl_d(pp%nrmax,0:pp%lmax)
  real(8) :: uvr(0:pp%lmax),duvr(0:pp%lmax)
  real(8) :: r,x,y,z
  real(8) :: xx
  real(8) :: ylm,dylm
  real(8) :: rshift(3)
  real(8) :: hvol

  hvol=hx*hy*hz

  if(iperiodic==0)then
    if(mod(lx(nl)-lx(1)+1,2)==1)then
      rshift(1)=0.d0
    else
      rshift(1)=-0.5d0*Hx
    end if
    if(mod(ly(nl)-ly(1)+1,2)==1)then
      rshift(2)=0.d0
    else
      rshift(2)=-0.5d0*Hy
    end if
    if(mod(lz(nl)-lz(1)+1,2)==1)then
      rshift(3)=0.d0
    else
      rshift(3)=-0.5d0*Hz
    end if
  else if(iperiodic==3)then 
    if(domain_parallel=='y')then
      rshift(1)=-Hx
      rshift(2)=-Hy
      rshift(3)=-Hz
    else
      rshift(:)=0.d0
    end if
  end if

  if(property /= 'update_wo_realloc') then

    do a=1,natom
      ik=kion(a)
      allocate(xn(0:pp%nrps(ik)-1),yn(0:pp%nrps(ik)-1),an(0:pp%nrps(ik)-2) &
              ,bn(0:pp%nrps(ik)-2),cn(0:pp%nrps(ik)-2),dn(0:pp%nrps(ik)-2))
  
      xn(0:pp%nrps(ik)-1) = pp%radnl(1:pp%nrps(ik),ik)
      do l=0,pp%mlps(ik)
        yn(0:pp%nrps(ik)-1) = pp%udvtbl(1:pp%nrps(ik),l,ik)
        call spline(pp%nrps(ik),xn,yn,an,bn,cn,dn)
        save_udvtbl_a(1:pp%nrps(ik)-1,l,a) = an(0:pp%nrps(ik)-2)
        save_udvtbl_b(1:pp%nrps(ik)-1,l,a) = bn(0:pp%nrps(ik)-2)
        save_udvtbl_c(1:pp%nrps(ik)-1,l,a) = cn(0:pp%nrps(ik)-2)
        save_udvtbl_d(1:pp%nrps(ik)-1,l,a) = dn(0:pp%nrps(ik)-2)
      end do
      deallocate(xn,yn,an,bn,cn,dn)
    enddo
  end if
  
  do a=1,natom
    ik=kion(a)

    if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
      allocate(xn(0:pp%nrps(ik)-1),yn(0:pp%nrps(ik)-1),an(0:pp%nrps(ik)-2) &
              ,bn(0:pp%nrps(ik)-2),cn(0:pp%nrps(ik)-2),dn(0:pp%nrps(ik)-2))
      xn(0:pp%nrps(ik)-1) = pp%radnl(1:pp%nrps(ik),ik)
      do l=0,pp%mlps(ik)
        yn(0:pp%nrps(ik)-1) = pp%dudvtbl(1:pp%nrps(ik),l,ik)
        call spline(pp%nrps(ik),xn,yn,an,bn,cn,dn)
        dudvtbl_a(1:pp%nrps(ik)-1,l) = an(0:pp%nrps(ik)-2)
        dudvtbl_b(1:pp%nrps(ik)-1,l) = bn(0:pp%nrps(ik)-2)
        dudvtbl_c(1:pp%nrps(ik)-1,l) = cn(0:pp%nrps(ik)-2)
        dudvtbl_d(1:pp%nrps(ik)-1,l) = dn(0:pp%nrps(ik)-2)
      end do
      deallocate(xn,yn,an,bn,cn,dn)
    endif
 
  !!$omp parallel
  !!$omp do private(j,x,y,z,r,ir,intr,xx,l,lm,m,uvr,duvr,ilma)
     do j=1,ppg%mps(a)
       x=ppg%jxyz(1,j,a)*hx+rshift(1)-(rion(1,a)+ppg%jxx(j,a)*alx)
       y=ppg%jxyz(2,j,a)*hy+rshift(2)-(rion(2,a)+ppg%jyy(j,a)*aly)
       z=ppg%jxyz(3,j,a)*hz+rshift(3)-(rion(3,a)+ppg%jzz(j,a)*alz)
       r=sqrt(x*x+y*y+z*z)+1d-50
       do ir=1,pp%nrps(ik)
         if(pp%radnl(ir,ik).gt.r) exit
       enddo
       intr=ir-1
       if (intr.lt.0.or.intr.ge.pp%nrps(ik))stop 'bad intr at prep_ps'
       xx = r - pp%radnl(intr,ik)
       do l=0,pp%mlps(ik)
         uvr(l) = save_udvtbl_a(intr,l,a)*xx**3 + save_udvtbl_b(intr,l,a)*xx**2 &
                + save_udvtbl_c(intr,l,a)*xx    + save_udvtbl_d(intr,l,a)
         if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
           duvr(l)=dudvtbl_a(intr,l)*xx**3 +dudvtbl_b(intr,l)*xx**2 &
                  +dudvtbl_c(intr,l)*xx    +dudvtbl_d(intr,l)
         endif
       enddo
 
       lm=0
       do l=0,pp%mlps(ik)
         if(pp%inorm(l,ik)==0) cycle
         do m=-l,l
           lm=lm+1
           ilma=ppg%lma_tbl(lm,a)
           ppg%uv(j,ilma)=uvr(l)*ylm(x,y,z,l,m)
           if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
             if(r>1d-6)then
               ppg%duv(j,ilma,1) = duvr(l)*(x/r)*ylm(x,y,z,l,m)+uvr(l)*dylm(x,y,z,l,m,1)
               ppg%duv(j,ilma,2) = duvr(l)*(y/r)*ylm(x,y,z,l,m)+uvr(l)*dylm(x,y,z,l,m,2)
               ppg%duv(j,ilma,3) = duvr(l)*(z/r)*ylm(x,y,z,l,m)+uvr(l)*dylm(x,y,z,l,m,3)
             else
               ppg%duv(j,ilma,1) = uvr(l)*dylm(x,y,z,l,m,1)
               ppg%duv(j,ilma,2) = uvr(l)*dylm(x,y,z,l,m,2)
               ppg%duv(j,ilma,3) = uvr(l)*dylm(x,y,z,l,m,3)
             end if
           end if
         enddo
       enddo
 
     enddo
 !!$omp end do
 !!$omp end parallel
 
  enddo

  lma=0
  do a=1,natom
    ik=kion(a)
    if(property /= 'update_wo_realloc') then
    do l=0,pp%mlps(ik)
      if(pp%inorm(l,ik)==0) cycle
      do m=-l,l
        lma=lma+1
        ppg%rinv_uvu(lma)=dble(pp%inorm(l,ik))*hvol
      enddo
    enddo
    endif
  enddo
 
end subroutine calc_uv

subroutine spline(Np,xn,yn,an,bn,cn,dn)
  integer,intent(in) :: Np
  real(8),intent(in) :: xn(0:Np-1),yn(0:Np-1)
  real(8),intent(out) :: an(0:Np-2),bn(0:Np-2),cn(0:Np-2),dn(0:Np-2)
  integer :: i,Npm2,info
  real(8) :: dxn(0:Np-1),dyn(0:Np-1),u(1:Np-2),v(1:Np-2),Amat(1:Np-2,1:Np-2)
  real(8) :: Amat_t(1:Np-2,1:Np-2)
! for lapack
  integer :: LWORK
  integer, allocatable :: IPIV(:) ! dimension N
  real(8), allocatable :: WORK(:) ! dimension LWORK
! for check inverse matrix problem
!  integer :: j,k
!  real(8) :: Amat_chk(1:Np-2,1:Np-2)
!  real(8) :: ss

  Npm2 = Np-2
  LWORK = Npm2*Npm2*6
  allocate(IPIV(Npm2),WORK(LWORK))


  do i = 0,Np-2
    dxn(i) = xn(i+1) - xn(i)
    dyn(i) = yn(i+1) - yn(i)
  end do

  do i = 1,Npm2
    v(i) = 6d0*(dyn(i)/dxn(i) - dyn(i-1)/dxn(i-1))
  end do

  Amat = 0d0
  Amat(1,1) = 2d0*(dxn(1) + dxn(0))
  Amat(1,2) = dxn(1)
  do i = 2,Npm2-1
    Amat(i,i+1) = dxn(i)
    Amat(i,i  ) = 2d0*(dxn(i)+dxn(i-1))
    Amat(i,i-1) = dxn(i-1)
  end do
  Amat(Npm2,Npm2  ) = 2d0*(dxn(Npm2)+dxn(Npm2-1))
  Amat(Npm2,Npm2-1) = dxn(Npm2-1)

! inverse matrix problem
  Amat_t = Amat


  call DGETRF(Npm2, Npm2, Amat_t, Npm2, IPIV, info)  ! factorize
  call DGETRI(Npm2, Amat_t, Npm2, IPIV, WORK, LWORK, info)  ! inverse

!  check inverse matrix problem
!  do i = 1,Npm2
!    do j = 1,Npm2
!      ss = 0d0
!      do k = 1,Npm2
!        ss = ss + Amat(i,k)*Amat_t(k,j)
!      end do
!      Amat_chk(i,j) = ss
!    end do
!  end do
!
!  do i = 1,Npm2
!    write(*,'(999e16.6e3)')(Amat_chk(i,j),j=1,Npm2)
!  end do
!
!  stop


  do i = 1,Npm2
    u(i) = sum(Amat_t(i,:)*v(:))
  end do

! for b
  bn(0) = 0d0
  bn(1:Np-2) = 0.5d0*u(1:Np-2)
! for a
  do i = 0,Npm2-1
    an(i) = (u(i+1) -2d0*bn(i))/(6d0*dxn(i))
  end do
  an(Npm2) = (0d0 -2d0*bn(Npm2))/(6d0*dxn(Npm2))
! for d
  dn(0:Npm2) = yn(0:Npm2)
! for c
  i=0
  cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*0.d0)/6d0
  do i = 1,Npm2-1
     cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*u(i))/6d0
  end do
  cn(Npm2) = dyn(Npm2)/dxn(Npm2) - dxn(Npm2)*(0d0+2d0*u(Npm2))/6d0

  return
end subroutine spline
