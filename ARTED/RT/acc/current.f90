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

! ####################
! current with OpenACC
! ####################

#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

#ifdef ARTED_CURRENT_PREPROCESSING
subroutine current_RT_preconditioner_LBLK(A, ikb_s,ikb_e)
  use Global_Variables
  use opt_variables
  implicit none
  complex(8),intent(in) :: A(0:NL-1, NBoccmax, NK_s:NK_e)
  integer,intent(in)    :: ikb_s,ikb_e

  call current_stencil_LBLK(A(:,:,:), ikb_s,ikb_e)
end subroutine
#endif

subroutine current_acc_impl(zutmp,jxs,jys,jzs)
  use Global_Variables
  use opt_variables
  implicit none
  complex(8),intent(in) :: zutmp(0:NL-1,NBoccmax,NK_s:NK_e)
  real(8),intent(out)   :: jxs,jys,jzs

  integer :: ikb,ib,ik,i,j,ix,iy,iz,ia
  real(8) :: kr,IaLxyz
  real(8) :: nabt(12)
  real(8) :: jx(NKB), jy(NKB), jz(NKB)
  integer :: ikb0, num_ikb1, ikb_s,ikb_e

  nabt( 1: 4) = nabx(1:4)
  nabt( 5: 8) = naby(1:4)
  nabt( 9:12) = nabz(1:4)

  IaLxyz = 1.0 / aLxyz

!$acc data pcopyin(zutmp) create(jx,jy,jz) pcopyout(ekr_omp) copyin(nabt) &
!$acc& pcopyin(jxyz,jxx,jyy,jzz,kAc,lx,ly,lz,Mps) 

!Constructing nonlocal part
!$acc kernels
!$acc loop collapse(2) independent gang
  do ik=NK_s,NK_e
  do ia=1,NI
!$acc loop independent vector(128)
  do j=1,Mps(ia)
    i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
    kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
    ekr_omp(j,ia,ik)=exp(zI*kr)
  end do
  end do
  end do
!$acc end kernels

  do ikb0 = 1, NKB, blk_nkb_current
    num_ikb1 = min(blk_nkb_current, NKB-ikb0+1)
    ikb_s = ikb0
    ikb_e = ikb0 + num_ikb1-1

    call init(zutmp(:,:,:), jx(:),jy(:),jz(:), ikb_s,ikb_e)
#ifndef ARTED_CURRENT_PREPROCESSING
    call current_stencil_LBLK(zutmp(:,:,:), ikb_s,ikb_e)
#endif
    call stencil(jx(:),jy(:),jz(:), ikb_s,ikb_e)
    call pseudo_pt(zutmp(:,:,:),IaLxyz,jx(:),jy(:),jz(:), ikb_s,ikb_e)
  end do

!$acc kernels
  jxs=0.d0
  jys=0.d0
  jzs=0.d0
!$acc loop gang vector reduction(+:jxs,jys,jzs)
  do ikb = 1, NKB
    jxs=jxs+jx(ikb)
    jys=jys+jy(ikb)
    jzs=jzs+jz(ikb)
  end do
!$acc end kernels

!$acc end data

contains
  subroutine init(zutmp,jx,jy,jz, ikb_s,ikb_e)
    use Global_Variables
    implicit none
    complex(8),intent(in) :: zutmp(0:NL-1, NBoccmax,NK_s:NK_e)
    real(8),intent(out)   :: jx(NKB),jy(NKB),jz(NKB)
    integer,intent(in)    :: ikb_s,ikb_e

    integer :: i
    real(8) :: jt
    integer :: ikb, ik,ib
!$acc data pcopy(jx,jy,jz) pcopyin(zutmp,occ,kAc, ik_table,ib_table)
!$acc kernels
!$acc loop gang
    do ikb = ikb_s, ikb_e
      ik=ik_table(ikb)
      ib=ib_table(ikb)

      jt=0.d0
!$acc loop vector(256) reduction(+:jt)
      do i=0,NL-1
        jt=jt+real(zutmp(i,ib,ik))**2+imag(zutmp(i,ib,ik))**2
      end do

      jx(ikb)=occ(ib,ik)*kAc(ik,1)*jt
      jy(ikb)=occ(ib,ik)*kAc(ik,2)*jt
      jz(ikb)=occ(ib,ik)*kAc(ik,3)*jt
    end do
!$acc end kernels
!$acc end data
  end subroutine

  subroutine stencil(jx,jy,jz, ikb_s,ikb_e)
    use Global_Variables
    use opt_variables
    implicit none
    real(8),intent(inout) :: jx(NKB),jy(NKB),jz(NKB)
    integer,intent(in)    :: ikb_s,ikb_e

    integer    :: ikb, ik,ib
!$acc data pcopy(jx,jy,jz) pcopyin(zcx,zcy,zcz,occ, ik_table,ib_table)
!$acc kernels
    do ikb = ikb_s, ikb_e
      ik=ik_table(ikb)
      ib=ib_table(ikb)

      jx(ikb)=jx(ikb)+zcx(ib,ik)*occ(ib,ik)
      jy(ikb)=jy(ikb)+zcy(ib,ik)*occ(ib,ik)
      jz(ikb)=jz(ikb)+zcz(ib,ik)*occ(ib,ik)
    end do
!$acc end kernels
!$acc end data
  end subroutine

  subroutine pseudo_pt(zutmp,IaLxyz,jx,jy,jz, ikb_s,ikb_e)
    use Global_Variables
    use omp_lib
    use opt_variables
    use timer
    implicit none
    complex(8),intent(in) :: zutmp(NL, NBoccmax,NK_s:NK_e)
    real(8),intent(in)    :: IaLxyz
    real(8),intent(out)   :: jx(NKB),jy(NKB),jz(NKB)
    integer,intent(in)    :: ikb_s,ikb_e

    integer    :: ilma,ia,j,i,ix,iy,iz
    real(8)    :: x,y,z
    real(8)    :: jxt,jyt,jzt
    complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz
    integer    :: ikb, ik,ib

!$acc data pcopy(jx,jy,jz) pcopyin(zutmp) create(t4cp_uVpsix,t4cp_uVpsiy,t4cp_uVpsiz) &
!$acc& pcopyin(ik_table,ib_table,a_tbl,Mps,Jxyz,Jxx,Jyy,Jzz,lx,ly,lz,uV,ekr_omp,iuV,occ)
!$acc kernels
!$acc loop collapse(2) gang
    do ikb = ikb_s, ikb_e
    do ilma=1,Nlma
      ik=ik_table(ikb)
      ib=ib_table(ikb)
      ia=a_tbl(ilma)
      uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
!$acc loop gang vector(256) reduction(+:uVpsi,uVpsix,uVpsiy,uVpsiz)
      do j=1,Mps(ia)
        i=Jxyz(j,ia)

        ix=Jxx(j,ia); x=Lx(i)*Hx-ix*aLx
        iy=Jyy(j,ia); y=Ly(i)*Hy-iy*aLy
        iz=Jzz(j,ia); z=Lz(i)*Hz-iz*aLz

        uVpsi =uVpsi +uV(j,ilma)*ekr_omp(j,ia,ik)  *zutmp(i,ib,ik)
        uVpsix=uVpsix+uV(j,ilma)*ekr_omp(j,ia,ik)*x*zutmp(i,ib,ik)
        uVpsiy=uVpsiy+uV(j,ilma)*ekr_omp(j,ia,ik)*y*zutmp(i,ib,ik)
        uVpsiz=uVpsiz+uV(j,ilma)*ekr_omp(j,ia,ik)*z*zutmp(i,ib,ik)
      end do
      uVpsi =uVpsi *Hxyz*iuV(ilma)
      uVpsix=uVpsix*Hxyz
      uVpsiy=uVpsiy*Hxyz
      uVpsiz=uVpsiz*Hxyz

      t4cp_uVpsix(ilma,ikb)=imag(conjg(uVpsix)*uVpsi)
      t4cp_uVpsiy(ilma,ikb)=imag(conjg(uVpsiy)*uVpsi)
      t4cp_uVpsiz(ilma,ikb)=imag(conjg(uVpsiz)*uVpsi)
    end do
    end do

!$acc loop gang
    do ikb = ikb_s, ikb_e
      ik=ik_table(ikb)
      ib=ib_table(ikb)

      jxt=0d0
      jyt=0d0
      jzt=0d0
!$acc loop vector(256) reduction(+:jxt,jyt,jzt)
      do ilma=1,Nlma
        jxt=jxt + t4cp_uVpsix(ilma,ikb)
        jyt=jyt + t4cp_uVpsiy(ilma,ikb)
        jzt=jzt + t4cp_uVpsiz(ilma,ikb)
      end do
      jxt=jxt * occ(ib,ik)*IaLxyz*2
      jyt=jyt * occ(ib,ik)*IaLxyz*2
      jzt=jzt * occ(ib,ik)*IaLxyz*2

      jx(ikb)=jx(ikb)*Hxyz*IaLxyz + jxt
      jy(ikb)=jy(ikb)*Hxyz*IaLxyz + jyt
      jz(ikb)=jz(ikb)*Hxyz*IaLxyz + jzt
    end do
!$acc end kernels
!$acc end data
  end subroutine
end subroutine
