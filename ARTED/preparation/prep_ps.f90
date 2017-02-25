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
!This file is "prep_ps.f90"
!This file conatain one soubroutine.
!SUBROUTINE prep_ps_periodic(property)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine prep_ps_periodic(property)
  use Global_Variables
  use communication
  implicit none
  character(11) :: property
  integer :: ik,n,i,a,j,ix,iy,iz,lma,l,m,lm,ir,intr
  integer :: lma_tbl((Lmax+1)**2,NI)
  real(8) :: G2sq,s,Vpsl_l(NL),G2,Gd,Gr,x,y,z,r
  real(8) :: Ylm,dYlm,uVr(0:Lmax),duVr(0:Lmax)
  complex(8) :: Vion_G(NG_s:NG_e)
  !spline interpolation
  real(8) :: xx
  real(8) :: udVtbl_a(Nrmax,0:Lmax),dudVtbl_a(Nrmax,0:Lmax)
  real(8) :: udVtbl_b(Nrmax,0:Lmax),dudVtbl_b(Nrmax,0:Lmax)
  real(8) :: udVtbl_c(Nrmax,0:Lmax),dudVtbl_c(Nrmax,0:Lmax)
  real(8) :: udVtbl_d(Nrmax,0:Lmax),dudVtbl_d(Nrmax,0:Lmax)
  real(8),allocatable :: xn(:),yn(:),an(:),bn(:),cn(:),dn(:)  
  real(8) :: vloc_av
  real(8) :: ratio1,ratio2,rc
  if(property == 'not_initial') then
    deallocate(a_tbl,uV,duV,iuV,Jxyz,Jxx,Jyy,Jzz)
    deallocate(ekr) ! sato
    deallocate(ekr_omp)
  end if

! local potential
!$omp parallel
!$omp do private(ik,n,G2sq,s,r,i,vloc_av) collapse(2)
  do ik=1,NE
    do n=NG_s,NG_e
      G2sq=sqrt(Gx(n)**2+Gy(n)**2+Gz(n)**2)
      s=0.d0
      if (n == nGzero) then
        do i=2,NRloc(ik)
           r=0.5d0*(rad(i,ik)+rad(i-1,ik))
           vloc_av = 0.5d0*(vloctbl(i,ik)+vloctbl(i-1,ik))
           s=s+4*Pi*(r**2*vloc_av+r*Zps(ik))*(rad(i,ik)-rad(i-1,ik))
        enddo
      else
        do i=2,NRloc(ik)
          r=0.5d0*(rad(i,ik)+rad(i-1,ik))
          vloc_av = 0.5d0*(vloctbl(i,ik)+vloctbl(i-1,ik))
          s=s+4*Pi*sin(G2sq*r)/(G2sq)*(r*vloc_av+Zps(ik))*(rad(i,ik)-rad(i-1,ik))
        enddo
      endif
      dVloc_G(n,ik)=s
    enddo
  enddo
!$omp end do
!$omp end parallel

  rhoion_G=0.d0
!$omp parallel private(a)
  do a=1,NI
!$omp do private(n)
    do n=NG_s,NG_e
      rhoion_G(n)=rhoion_G(n)+Zps(Kion(a))/aLxyz*exp(-zI*(Gx(n)*Rion(1,a)+Gy(n)*Rion(2,a)+Gz(n)*Rion(3,a)))
    enddo
!$omp end do
  enddo
!$omp end parallel

  Vion_G=0.d0
!$omp parallel private(a,ik)
  do a=1,NI
    ik=Kion(a)
!$omp do private(n,G2,Gd)
    do n=NG_s,NG_e
      G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
      Gd=Gx(n)*Rion(1,a)+Gy(n)*Rion(2,a)+Gz(n)*Rion(3,a)
      Vion_G(n)=Vion_G(n)+dVloc_G(n,ik)*exp(-zI*Gd)/aLxyz
      if(n == nGzero) cycle
      Vion_G(n)=Vion_G(n)-4*Pi/G2*Zps(ik)*exp(-zI*Gd)/aLxyz
    enddo
!$omp end do
  enddo
!$omp end parallel

  Vpsl_l=0.d0
!$omp parallel private(n)
  do n=NG_s,NG_e
!$omp do private(i,Gr)
    do i=1,NL
      Gr=Gx(n)*Lx(i)*Hx+Gy(n)*Ly(i)*Hy+Gz(n)*Lz(i)*Hz
      Vpsl_l(i)=Vpsl_l(i)+Vion_G(n)*exp(zI*Gr)
    enddo
!$omp end do
  enddo
!$omp end parallel

  call comm_summation(Vpsl_l,Vpsl,NL,proc_group(2))

! nonlocal potential
  if (comm_is_root() .and. property=='initial') then
    write(*,*) ''
    write(*,*) '============nonlocal grid data=============='
  endif
  do a=1,NI
    ik=Kion(a)
    j=0
    do ix=-2,2
    do iy=-2,2
    do iz=-2,2
      do i=1,NL
        x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
        y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
        z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
        r=sqrt(x**2+y**2+z**2)
        if (r<Rps(ik)) then
          j=j+1
        endif
      enddo
    enddo
    enddo
    enddo
    Mps(a)=j
    if (comm_is_root() .and. property == 'initial') then
      write(*,*) 'a =',a,'Mps(a) =',Mps(a)
    endif
  end do
  Nps=maxval(Mps(:))
  allocate(Jxyz(Nps,NI),Jxx(Nps,NI),Jyy(Nps,NI),Jzz(Nps,NI))
  allocate(ekr(Nps,NI)) ! sato
  allocate(ekr_omp(Nps,NI,NK_s:NK_e))

  do a=1,NI
    ik=Kion(a)
    j=0
    do ix=-2,2
    do iy=-2,2
    do iz=-2,2
      do i=1,NL
        x=Lx(i)*Hx-(Rion(1,a)+ix*aLx)
        y=Ly(i)*Hy-(Rion(2,a)+iy*aLy)
        z=Lz(i)*Hz-(Rion(3,a)+iz*aLz)
        r=sqrt(x**2+y**2+z**2)
        if (r<Rps(ik)) then
          j=j+1
          if (j<=Nps) then
            Jxyz(j,a)=i
            Jxx(j,a)=ix
            Jyy(j,a)=iy
            Jzz(j,a)=iz
          endif
        endif
      enddo
    enddo
    enddo
    enddo
  end do

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

  allocate(a_tbl(Nlma),uV(Nps,Nlma),iuV(Nlma),duV(Nps,Nlma,3))

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

  do a=1,NI
    ik=Kion(a)

    allocate(xn(0:NRps(ik)-1),yn(0:NRps(ik)-1),an(0:NRps(ik)-2) &
         ,bn(0:NRps(ik)-2),cn(0:NRps(ik)-2),dn(0:NRps(ik)-2))
    
    xn(0:NRps(ik)-1) = radnl(1:NRps(ik),ik)
    do l=0,Mlps(ik)
       yn(0:NRps(ik)-1) = udVtbl(1:NRps(ik),l,ik)
       call spline(NRps(ik),xn,yn,an,bn,cn,dn)
       udVtbl_a(1:NRps(ik)-1,l) = an(0:NRps(ik)-2)
       udVtbl_b(1:NRps(ik)-1,l) = bn(0:NRps(ik)-2)
       udVtbl_c(1:NRps(ik)-1,l) = cn(0:NRps(ik)-2)
       udVtbl_d(1:NRps(ik)-1,l) = dn(0:NRps(ik)-2)

       yn(0:NRps(ik)-1) = dudVtbl(1:NRps(ik),l,ik)
       call spline(NRps(ik),xn,yn,an,bn,cn,dn)
       dudVtbl_a(1:NRps(ik)-1,l) = an(0:NRps(ik)-2)
       dudVtbl_b(1:NRps(ik)-1,l) = bn(0:NRps(ik)-2)
       dudVtbl_c(1:NRps(ik)-1,l) = cn(0:NRps(ik)-2)
       dudVtbl_d(1:NRps(ik)-1,l) = dn(0:NRps(ik)-2)        
    end do
    
    do j=1,Mps(a)
      x=Lx(Jxyz(j,a))*Hx-(Rion(1,a)+Jxx(j,a)*aLx)
      y=Ly(Jxyz(j,a))*Hy-(Rion(2,a)+Jyy(j,a)*aLy)
      z=Lz(Jxyz(j,a))*Hz-(Rion(3,a)+Jzz(j,a)*aLz)
      r=sqrt(x**2+y**2+z**2)+1d-50
      do ir=1,NRps(ik)
        if(radnl(ir,ik).gt.r) exit
      enddo
      intr=ir-1
      if (intr.lt.0.or.intr.ge.NRps(ik))stop 'bad intr at prep_ps'
      xx = r - radnl(intr,ik) 
      do l=0,Mlps(ik)
         uVr(l) = udVtbl_a(intr,l)*xx**3 + udVtbl_b(intr,l)*xx**2 &
                 +udVtbl_c(intr,l)*xx    + udVtbl_d(intr,l)
         duVr(l)=dudVtbl_a(intr,l)*xx**3 +dudVtbl_b(intr,l)*xx**2 &
                +dudVtbl_c(intr,l)*xx    +dudVtbl_d(intr,l)         
      enddo
      lm=0
      do l=0,Mlps(ik)
        if(inorm(l,ik)==0) cycle
        do m=-l,l
          lm=lm+1
          uV(j,lma_tbl(lm,a))=uVr(l)*Ylm(x,y,z,l,m)
          if(r>1d-6)then
             duV(j,lma_tbl(lm,a),1)=duVr(l)*(x/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,1)
             duV(j,lma_tbl(lm,a),2)=duVr(l)*(y/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,2)
             duV(j,lma_tbl(lm,a),3)=duVr(l)*(z/r)*Ylm(x,y,z,l,m)+uVr(l)*dYlm(x,y,z,l,m,3)
          else
             duV(j,lma_tbl(lm,a),1)=uVr(l)*dYlm(x,y,z,l,m,1)
             duV(j,lma_tbl(lm,a),2)=uVr(l)*dYlm(x,y,z,l,m,2)
             duV(j,lma_tbl(lm,a),3)=uVr(l)*dYlm(x,y,z,l,m,3)
          end if
        enddo
      enddo
    enddo
    lm=0
    do l=0,Mlps(ik)
      if(inorm(l,ik)==0) cycle
      do m=-l,l
        lm=lm+1
        iuV(lma_tbl(lm,a))=inorm(l,ik)
      enddo
    enddo
    deallocate(xn,yn,an,bn,cn,dn)
  enddo

! nonlinear core-correction
  rho_nlcc = 0d0
  tau_nlcc = 0d0
  if(flag_nlcc)then
    if(comm_is_root())write(*,"(A)")"Preparation: Non-linear core correction"
    do a=1,NI
      ik=Kion(a)
      rc = 15d0 ! maximum
      do i=1,Nrmax
        if(rho_nlcc_tbl(i,ik) + tau_nlcc_tbl(i,ik) < 1d-6)then
          rc = rad(i,ik)
          exit
        end if
        if(i == Nrmax) stop"no-cut-off"
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

  return
End Subroutine prep_ps_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
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
    Amat(i,i) = 2d0*(dxn(i)+dxn(i-1))
    Amat(i,i-1) = dxn(i-1)
  end do
  Amat(Npm2,Npm2) = 2d0*(dxn(Npm2)+dxn(Npm2-1))
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
!    write(*,*)u(i),v(i)
  end do
!  stop


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
  do i = 0,Npm2-1
    if (i==0) then
      cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*0.d0)/6d0
    else
      cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*u(i))/6d0 ! OOR when i=0
    endif
  end do
  cn(Npm2) = dyn(Npm2)/dxn(Npm2) - dxn(Npm2)*(0d0+2d0*u(Npm2))/6d0

  return
end subroutine spline
