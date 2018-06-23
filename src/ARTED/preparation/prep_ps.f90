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
!This file is "prep_ps.f90"
!This file conatain one soubroutine.
!SUBROUTINE prep_ps_periodic(property)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine prep_ps_periodic(property)
  use Global_Variables
  use salmon_parallel, only: nproc_id_global, nproc_group_tdks
  use salmon_communication, only: comm_summation, comm_is_root
  use salmon_math
  use opt_variables, only: zJxyz,zKxyz,init_for_padding, nprojector,idx_proj,idx_lma,pseudo_start_idx, init_projector
#ifdef ARTED_LBLK
  use opt_variables, only: t4ppt_nlma,t4ppt_nlma,t4ppt_i2vi,t4ppt_vi2i,t4ppt_ilma,t4ppt_j,  opt_vars_init_t4ppt
#endif
  implicit none
  character(17) :: property
  logical :: flag_alloc1, flag_alloc2
  integer :: ik,n,i,a,j,ix,iy,iz,lma,l,m,lm,ir,intr
  integer :: PNLx,PNLy,PNLz,ilma,narray
  real(8) :: G2sq,s,G2,Gd,Gr,x,y,z,r,dr,tmpx,tmpy,tmpz
  real(8) :: Ylm,dYlm,uVr(0:Lmax),duVr(0:Lmax),Vpsl_ia_l(NL,NI) !,Vpsl_l(NL)
  complex(8) :: Vion_G_ia(NG_s:NG_e,NI),tmp_exp !, Vion_G(NG_s:NG_e)
  !spline interpolation
  real(8) :: xx
  real(8) :: dudVtbl_a(Nrmax,0:Lmax), dudVtbl_b(Nrmax,0:Lmax)
  real(8) :: dudVtbl_c(Nrmax,0:Lmax), dudVtbl_d(Nrmax,0:Lmax)
!  real(8) :: udVtbl_a(Nrmax,0:Lmax), udVtbl_b(Nrmax,0:Lmax)
!  real(8) :: udVtbl_c(Nrmax,0:Lmax), udVtbl_d(Nrmax,0:Lmax)
  real(8),allocatable :: xn(:),yn(:),an(:),bn(:),cn(:),dn(:)  
  real(8) :: vloc_av, ratio1,ratio2,rc


  !(Local pseudopotential in G-space (radial part))
  if(property == 'initial') then

    allocate(lma_tbl((Lmax+1)**2,NI))

!$omp parallel
!$omp do private(ik,n,G2sq,s,r,dr,i,vloc_av) collapse(2)
    do ik=1,NE
       do n=NG_s,NG_e
          G2sq=sqrt(Gx(n)**2+Gy(n)**2+Gz(n)**2)
          s=0.d0
          if (n == nGzero) then
             do i=2,NRloc(ik)
                 r=0.5d0*(rad(i,ik)+rad(i-1,ik))
                dr=rad(i,ik)-rad(i-1,ik)
                vloc_av = 0.5d0*(vloctbl(i,ik)+vloctbl(i-1,ik))
                s=s+4d0*Pi*(r**2*vloc_av+r*Zps(ik))*dr
             enddo
          else
             do i=2,NRloc(ik)
                 r=0.5d0*(rad(i,ik)+rad(i-1,ik))
                dr=rad(i,ik)-rad(i-1,ik)
                vloc_av = 0.5d0*(vloctbl(i,ik)+vloctbl(i-1,ik))
                s=s+4d0*Pi*sin(G2sq*r)/G2sq*(r*vloc_av+Zps(ik))*dr !Vloc - coulomb
             enddo
          endif
          dVloc_G(n,ik)=s
       enddo
    enddo
!$omp end do
!$omp end parallel

    !(this save is used only for opt and md options)
    save_dVloc_G(:,:)=dVloc_G(:,:)
     
  else

    dVloc_G(:,:)=save_dVloc_G(:,:)

    if(property == 'update_all') then
#define SAFE_DEALLOCATE(var) if(allocated(var)) deallocate(var)
#ifdef ARTED_LBLK
       SAFE_DEALLOCATE(t4ppt_nlma)
       SAFE_DEALLOCATE(t4ppt_i2vi)
       SAFE_DEALLOCATE(t4ppt_vi2i)
       SAFE_DEALLOCATE(t4ppt_ilma)
       SAFE_DEALLOCATE(t4ppt_j)
#endif
       SAFE_DEALLOCATE(nprojector)
       SAFE_DEALLOCATE(idx_proj)
       SAFE_DEALLOCATE(idx_lma)
       SAFE_DEALLOCATE(pseudo_start_idx)
    endif
  endif


  !(Local pseudopotential: Vlocal in G-space(=Vion_G))
  Vion_G_ia=0.d0
 !Vion_G   =0.d0
  rhoion_G =0.d0
!$omp parallel private(a,ik)
  do a=1,NI
    ik=Kion(a)
!$omp do private(n,G2,Gd,tmp_exp)
    do n=NG_s,NG_e
      Gd=Gx(n)*Rion(1,a)+Gy(n)*Rion(2,a)+Gz(n)*Rion(3,a)
      tmp_exp = exp(-zI*Gd)/aLxyz
     !Vion_G(n)     = Vion_G(n)      + dVloc_G(n,ik)*tmp_exp
      Vion_G_ia(n,a)= Vion_G_ia(n,a) + dVloc_G(n,ik)*tmp_exp
      rhoion_G(n)   = rhoion_G(n) + Zps(ik)*tmp_exp
      if(n == nGzero) cycle
      !(add coulomb as dVloc_G is given by Vloc - coulomb)
      G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
     !Vion_G(n)     = Vion_G(n)      -4d0*Pi/G2*Zps(ik)*tmp_exp
      Vion_G_ia(n,a)= Vion_G_ia(n,a) -4d0*Pi/G2*Zps(ik)*tmp_exp
    enddo
!$omp end do
  enddo
!$omp end parallel


  !(Local pseudopotential: Vlocal(=Vpsl) in real-space)
  Vpsl_ia_l=0.d0
 !Vpsl_l   =0.d0
!$omp parallel private(n)
  do n=NG_s,NG_e
!$omp do private(i,Gr,a,tmp_exp)
     do i=1,NL
        Gr = Gx(n)*Lx(i)*Hx+Gy(n)*Ly(i)*Hy+Gz(n)*Lz(i)*Hz
        tmp_exp = exp(zI*Gr)
       !Vpsl_l(i) = Vpsl_l(i) + Vion_G(n)*tmp_exp
        do a=1,NI
           Vpsl_ia_l(i,a)= Vpsl_ia_l(i,a)+ Vion_G_ia(n,a)*tmp_exp
        enddo
     enddo
!$omp end do
  enddo
!$omp end parallel

 !call comm_summation(Vpsl_l,Vpsl,NL,nproc_group_tdks)
  call comm_summation(Vpsl_ia_l,Vpsl_ia,NL*NI,nproc_group_tdks)

!$omp parallel
!$omp do private(i)
  do i=1,NL
     Vpsl(i) = sum(Vpsl_ia(i,:))
  enddo
!$omp end do
!$omp end parallel

  !(Non-Local pseudopotential)
  if(property /= 'update_wo_realloc') then

  if (comm_is_root(nproc_id_global) .and. property=='initial') then
    write(*,*) ''
    write(*,*) '============nonlocal grid data=============='
  endif

!$omp parallel
!$omp do private(a,ik,j,i,ix,iy,iz,x,y,z,r,tmpx,tmpy,tmpz)
  do a=1,NI
     ik=Kion(a)
     j=0
     do ix=-2,2
     do iy=-2,2
     do iz=-2,2
        tmpx = Rion(1,a)+ix*aLx
        tmpy = Rion(2,a)+iy*aLy
        tmpz = Rion(3,a)+iz*aLz
        do i=1,NL
           x=Lx(i)*Hx-tmpx
           y=Ly(i)*Hy-tmpy
           z=Lz(i)*Hz-tmpz
           r=sqrt(x*x+y*y+z*z)
           if (r<Rps(ik)) j=j+1
        enddo
     enddo
     enddo
     enddo
    Mps(a)=j
  end do
!$omp end do
!$omp end parallel

  Nps=maxval(Mps(:))

  if (comm_is_root(nproc_id_global) .and. property == 'initial') then
     do a=1,NI
        write(*,*) 'a =',a,'Mps(a) =',Mps(a)
     end do
  endif

  !(allocate/deallocate with Nps)
  if(property == 'initial') then
     flag_alloc1=.true.
  else if(property == 'update_all') then
     narray=ubound(Jxyz,1)
     if(Nps.ne.narray)then
        deallocate(Jxyz,Jxx,Jyy,Jzz,zJxyz)
        deallocate(ekr,ekr_omp)
#ifdef ARTED_STENCIL_PADDING
        deallocate(zKxyz)
#endif
        flag_alloc1=.true.
     else
        flag_alloc1=.false.
     endif
  endif
  if(flag_alloc1)then
     allocate(Jxyz(Nps,NI),Jxx(Nps,NI),Jyy(Nps,NI),Jzz(Nps,NI),zJxyz(Nps,NI))
     allocate(ekr_omp(Nps,NI,NK_s:NK_e),ekr(Nps,NI))
#ifdef ARTED_STENCIL_PADDING
     allocate(zKxyz(Nps,NI))
#endif
  endif

!$omp parallel
!$omp do private(a,ik,j,ix,iy,iz,tmpx,tmpy,tmpz,i,x,y,z,r)
  do a=1,NI
     ik=Kion(a)
     j=0
     do ix=-2,2
     do iy=-2,2
     do iz=-2,2
        tmpx = Rion(1,a)+ix*aLx
        tmpy = Rion(2,a)+iy*aLy
        tmpz = Rion(3,a)+iz*aLz
        do i=1,NL
           x=Lx(i)*Hx-tmpx
           y=Ly(i)*Hy-tmpy
           z=Lz(i)*Hz-tmpz
           r=sqrt(x*x+y*y+z*z)
           if (r<Rps(ik)) then
              j=j+1
              if (j<=Nps) then
                 Jxyz(j,a)=i
                 Jxx( j,a)=ix
                 Jyy( j,a)=iy
                 Jzz( j,a)=iz
              endif
           endif
        enddo
     enddo
     enddo
     enddo
  end do
!$omp end do
!$omp end parallel

  if(property == 'update_all') then
     zJxyz(1:Nps,1:NI) = Jxyz(1:Nps,1:NI) - 1

#ifdef ARTED_STENCIL_PADDING
     !call init_for_padding
     PNLx = NLx
     PNLy = NLy + 1
     PNLz = NLz

!$omp parallel
!$omp do private(a,j,i)
    do a=1,NI
    do j=1,Mps(a)
       i=Jxyz(j,a)
       zKxyz(j,a)=Lx(i)*PNLy*PNLz + Ly(i)*PNLz + Lz(i)
    enddo
    enddo
!$omp end do
!$omp end parallel
#endif

  endif

  lma=0
  do a=1,NI
    ik=Kion(a)
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lma=lma+1
      enddo
    enddo
  enddo
  Nlma=lma


  !(allocate/deallocate with Nlma)
  if(property == 'initial') then
     flag_alloc2=.true.
  else if(property == 'update_all') then
     narray=ubound(a_tbl,1)
     if(Nlma.ne.narray .or. flag_alloc1)then
        deallocate(a_tbl,uV,duV,iuV,zproj)
        flag_alloc2=.true.
     else
        flag_alloc2=.false.
     endif
  endif
  if(flag_alloc2)then
     allocate(a_tbl(Nlma),uV(Nps,Nlma),iuV(Nlma),duV(Nps,Nlma,3))
     allocate(zproj(Nps,Nlma,NK_s:NK_e))
  endif

  lma=0
  do a=1,NI
    ik=Kion(a)
    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        lma=lma+1
        a_tbl(lma)=a
        lma_tbl(lm,a)=lma
      enddo
    enddo
  enddo

  endif  !for /= 'update_wo_realloc'


  if(property /= 'update_wo_realloc') then

  if(property == 'initial') then
     allocate( save_udVtbl_a(Nrmax,0:Lmax,NI) )
     allocate( save_udVtbl_b(Nrmax,0:Lmax,NI) )
     allocate( save_udVtbl_c(Nrmax,0:Lmax,NI) )
     allocate( save_udVtbl_d(Nrmax,0:Lmax,NI) )
  endif

  narray=-99
  do a=1,NI
    ik=Kion(a)

    if(a.ne.1)narray=ubound(xn,1)
    if(narray.ne.NRps(ik)-1) then
       if(a.ne.1) deallocate(xn,yn,an,bn,cn,dn)
       allocate(xn(0:NRps(ik)-1),yn(0:NRps(ik)-1),an(0:NRps(ik)-2) &
               ,bn(0:NRps(ik)-2),cn(0:NRps(ik)-2),dn(0:NRps(ik)-2))
    endif

    xn(0:NRps(ik)-1) = radnl(1:NRps(ik),ik)
    do l=0,Mlps(ik)
       yn(0:NRps(ik)-1) = udVtbl(1:NRps(ik),l,ik)
       call spline(NRps(ik),xn,yn,an,bn,cn,dn)
       save_udVtbl_a(1:NRps(ik)-1,l,a) = an(0:NRps(ik)-2)
       save_udVtbl_b(1:NRps(ik)-1,l,a) = bn(0:NRps(ik)-2)
       save_udVtbl_c(1:NRps(ik)-1,l,a) = cn(0:NRps(ik)-2)
       save_udVtbl_d(1:NRps(ik)-1,l,a) = dn(0:NRps(ik)-2)
    end do

  enddo
  endif


  narray=-99
  do a=1,NI
    ik=Kion(a)

    if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
    if(a.ne.1)narray=ubound(xn,1)
    if(narray.ne.NRps(ik)-1) then
       if(a.ne.1) deallocate(xn,yn,an,bn,cn,dn)
       allocate(xn(0:NRps(ik)-1),yn(0:NRps(ik)-1),an(0:NRps(ik)-2) &
               ,bn(0:NRps(ik)-2),cn(0:NRps(ik)-2),dn(0:NRps(ik)-2))
    endif
    xn(0:NRps(ik)-1) = radnl(1:NRps(ik),ik)
    do l=0,Mlps(ik)
       !yn(0:NRps(ik)-1) = udVtbl(1:NRps(ik),l,ik)
       !call spline(NRps(ik),xn,yn,an,bn,cn,dn)
       !udVtbl_a(1:NRps(ik)-1,l) = an(0:NRps(ik)-2)
       !udVtbl_b(1:NRps(ik)-1,l) = bn(0:NRps(ik)-2)
       !udVtbl_c(1:NRps(ik)-1,l) = cn(0:NRps(ik)-2)
       !udVtbl_d(1:NRps(ik)-1,l) = dn(0:NRps(ik)-2)
       yn(0:NRps(ik)-1) = dudVtbl(1:NRps(ik),l,ik)
       call spline(NRps(ik),xn,yn,an,bn,cn,dn)
       dudVtbl_a(1:NRps(ik)-1,l) = an(0:NRps(ik)-2)
       dudVtbl_b(1:NRps(ik)-1,l) = bn(0:NRps(ik)-2)
       dudVtbl_c(1:NRps(ik)-1,l) = cn(0:NRps(ik)-2)
       dudVtbl_d(1:NRps(ik)-1,l) = dn(0:NRps(ik)-2)        
    end do
    endif

!$omp parallel
!$omp do private(j,x,y,z,r,ir,intr,xx,l,lm,m,uVr,duVr,ilma)
    do j=1,Mps(a)
      x=Lx(Jxyz(j,a))*Hx-(Rion(1,a)+Jxx(j,a)*aLx)
      y=Ly(Jxyz(j,a))*Hy-(Rion(2,a)+Jyy(j,a)*aLy)
      z=Lz(Jxyz(j,a))*Hz-(Rion(3,a)+Jzz(j,a)*aLz)
      r=sqrt(x*x+y*y+z*z)+1d-50
      do ir=1,NRps(ik)
        if(radnl(ir,ik).gt.r) exit
      enddo
      intr=ir-1
     !if (intr.lt.0.or.intr.ge.NRps(ik))stop 'bad intr at prep_ps'
      xx = r - radnl(intr,ik) 
      do l=0,Mlps(ik)
         uVr(l) = save_udVtbl_a(intr,l,a)*xx**3 + save_udVtbl_b(intr,l,a)*xx**2 &
                + save_udVtbl_c(intr,l,a)*xx    + save_udVtbl_d(intr,l,a)
        !uVr(l) = udVtbl_a(intr,l)*xx**3 + udVtbl_b(intr,l)*xx**2 &
        !        +udVtbl_c(intr,l)*xx    + udVtbl_d(intr,l)
         if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
         duVr(l)=dudVtbl_a(intr,l)*xx**3 +dudVtbl_b(intr,l)*xx**2 &
                +dudVtbl_c(intr,l)*xx    +dudVtbl_d(intr,l)         
         endif
      enddo

      lm=0
      do l=0,Mlps(ik)
        if(inorm(l,ik)==0) cycle
        do m=-l,l
          lm=lm+1
          ilma=lma_tbl(lm,a)
          uV(j,ilma)=uVr(l)*Ylm(x,y,z,l,m)
          if(.not.flag_use_grad_wf_on_force)then !legacy for ion-force
          if(r>1d-6)then
             duV(j,ilma,1) = duVr(l)*(x/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,1)
             duV(j,ilma,2) = duVr(l)*(y/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,2)
             duV(j,ilma,3) = duVr(l)*(z/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,3)
          else
             duV(j,ilma,1) = uVr(l)*dYlm(x,y,z,l,m,1)
             duV(j,ilma,2) = uVr(l)*dYlm(x,y,z,l,m,2)
             duV(j,ilma,3) = uVr(l)*dYlm(x,y,z,l,m,3)
          end if
          endif
        enddo
      enddo

    enddo
!$omp end do
!$omp end parallel

    if(property /= 'update_wo_realloc') then
    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        iuV(lma_tbl(lm,a))=inorm(l,ik)
      enddo
    enddo
    endif

  enddo

! nonlinear core-correction
  rho_nlcc = 0d0
  tau_nlcc = 0d0
  if(flag_nlcc)then
    if(comm_is_root(nproc_id_global))write(*,"(A)")"Preparation: Non-linear core correction"
    do a=1,NI
      ik=Kion(a)
      rc = 15d0 ! maximum
      do i=1,Nrmax
        if(rho_nlcc_tbl(i,ik) + tau_nlcc_tbl(i,ik) < 1d-6)then
          rc = rad(i,ik)
          exit
        end if
        if(i == Nrmax) stop "no-cut-off"
      end do
      
      do ix=-2,2; do iy=-2,2; do iz=-2,2
        do i=1,NL
          x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
          y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
          z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
          r=sqrt(x**2+y**2+z**2)
          if(r > rc)cycle
          
          do ir=1,NRmax
            if(rad(ir,ik).gt.r) exit
          enddo
          intr=ir-1
          if (intr.lt.0.or.intr.ge.NRmax)stop 'bad intr at prep_ps'
          ratio1=(r-rad(intr,ik))/(rad(intr+1,ik)-rad(intr,ik))
          ratio2=1-ratio1
          rho_nlcc(i) = rho_nlcc(i) &
            +ratio1*rho_nlcc_tbl(intr+1,ik)+ratio2*rho_nlcc_tbl(intr,ik)
          tau_nlcc(i) = tau_nlcc(i) &
            +ratio1*tau_nlcc_tbl(intr+1,ik)+ratio2*tau_nlcc_tbl(intr,ik)
          
        enddo
        
      end do; end do; end do
    end do
  end if


  if(property == 'update_all') then
#ifdef ARTED_STENCIL_PADDING
    call init_projector(zKxyz)
#else
    call init_projector(zJxyz)
#endif
#ifdef ARTED_LBLK
     call opt_vars_init_t4ppt
#endif
  endif

  return
End Subroutine prep_ps_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
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
