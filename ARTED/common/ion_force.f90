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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine Ion_Force_omp(Rion_update,GS_RT)
  use Global_Variables, only: zu,zu_GS,NB,NBoccmax
  implicit none
  character(2),intent(in) :: GS_RT
  character(3),intent(in) :: Rion_update

  if(GS_RT == 'GS') then
    call impl(Rion_update,zu_GS,NB)
  else if(GS_RT == 'RT') then
    call impl(Rion_update,zu,NBoccmax)
  end if

contains
  subroutine impl(Rion_update,zutmp,zu_NB)
    use Global_Variables
    use communication
    use timer
    implicit none
    character(3),intent(in)  :: Rion_update
    integer,intent(in)       :: zu_NB
    complex(8),intent(inout) :: zutmp(NL,zu_NB,NK_s:NK_e)

    integer      :: ia,ib,ilma,ik,ix,iy,iz,n,j,i
    real(8)      :: rab(3),rab2,Gvec(3),G2,Gd,ftmp_l(3,NI),kr
    complex(8)   :: uVpsi,duVpsi(3)
    real(8)      :: ftmp_l_kl(3,NI,NK_s:NK_e)

    call timer_begin(LOG_ION_FORCE)

    !ion
    if (Rion_update == 'on') then
      ftmp_l=0.d0
!$omp parallel
      do ia=1,NI
!$omp do private(ik,ix,iy,iz,ib,rab,rab2) collapse(4)
      do ix=-NEwald,NEwald
      do iy=-NEwald,NEwald
      do iz=-NEwald,NEwald
      do ib=1,NI
        ik=Kion(ia)
        if(ix**2+iy**2+iz**2 == 0 .and. ia == ib) cycle
        rab(1)=Rion(1,ia)-ix*aLx-Rion(1,ib)
        rab(2)=Rion(2,ia)-iy*aLy-Rion(2,ib)
        rab(3)=Rion(3,ia)-iz*aLz-Rion(3,ib)
        rab2=sum(rab(:)**2)
        ftmp_l(:,ia)=ftmp_l(:,ia)&
             &-Zps(Kion(ia))*Zps(Kion(ib))*rab(:)/sqrt(rab2)*(-erfc(sqrt(aEwald*rab2))/rab2&
             &-2*sqrt(aEwald/(rab2*Pi))*exp(-aEwald*rab2))
      end do
      end do
      end do
      end do
!$omp end do
      end do
!$omp end parallel
      Fion=ftmp_l
    end if

    ftmp_l=0.d0
    ftmp_l_kl=0.d0

!$omp parallel private(ia)

    !loc
    do ia=1,NI
!$omp do private(ik,n,Gvec,G2,Gd)
    do n=NG_s,NG_e
      if(n == nGzero) cycle
      ik=Kion(ia)
      Gvec(1)=Gx(n); Gvec(2)=Gy(n); Gvec(3)=Gz(n)
      G2=sum(Gvec(:)**2)
      Gd=sum(Gvec(:)*Rion(:,ia))
      ftmp_l(:,ia)=ftmp_l(:,ia)&
           &+zI*Gvec(:)*(4*Pi/G2)*Zps(ik)*(rhoe_G(n)*exp(zI*Gd) &
           &+0.5d0*exp(-G2/(4*aEwald))*(conjg(rhoion_G(n))*exp(-zI*Gd)-rhoion_G(n)*exp(zI*Gd)))&
           &+conjg(rhoe_G(n))*dVloc_G(n,ik)*zI*Gvec(:)*exp(-zI*Gd)
    end do
!$omp end do
    end do

    !nonlocal
!$omp do private(ik,j,i,ix,iy,iz,kr) collapse(2)
    do ik=NK_s,NK_e
    do ia=1,NI
    do j=1,Mps(ia)
      i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
      kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
      ekr_omp(j,ia,ik)=exp(zI*kr)
    end do
    end do
    end do
!$omp end do

!$omp do private(ik,j,i,ib,ilma,uVpsi,duVpsi) collapse(2)
    do ik=NK_s,NK_e
    do ib=1,NBoccmax
      do ilma=1,Nlma
        ia=a_tbl(ilma)
        uVpsi=0.d0
        duVpsi(:)=0.d0
        do j=1,Mps(ia)
          i=Jxyz(j,ia)
          uVpsi=uVpsi+uV(j,ilma)*ekr_omp(j,ia,ik)*zutmp(i,ib,ik)
          duVpsi(:)=duVpsi(:)+duV(j,ilma,:)*ekr_omp(j,ia,ik)*zutmp(i,ib,ik)
        end do
        uVpsi=uVpsi*Hxyz; duVpsi(:)=duVpsi(:)*Hxyz
        ftmp_l_kl(:,ia,ik)=ftmp_l_kl(:,ia,ik)+(conjg(uVpsi)*duVpsi(:)+uVpsi*conjg(duVpsi(:)))*iuV(ilma)*occ(ib,ik)
      end do
    end do
    end do
!$omp end do

!$omp end parallel

    ftmp_l(:,:)=ftmp_l(:,:)+ftmp_l_kl(:,:,NK_s)
    do ik=NK_s+1,NK_e
      ftmp_l(:,:)=ftmp_l(:,:)+ftmp_l_kl(:,:,ik)
    end do

    call timer_end(LOG_ION_FORCE)

    call timer_begin(LOG_ALLREDUCE)
    call comm_summation(ftmp_l,fnl,3*NI,proc_group(2))
    force=Floc+Fnl+Fion
    call timer_end(LOG_ALLREDUCE)
  end subroutine
end subroutine Ion_Force_omp
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
