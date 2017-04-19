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
!This file is "current.f90"
!This file contain one subroutine.
!SUBROUTINE current
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif


#ifdef ARTED_CURRENT_PREPROCESSING
subroutine current_RT_preconditioner(ib,ik,A)
  use Global_Variables
  use opt_variables
  implicit none
  integer,intent(in)    :: ib,ik
  complex(8),intent(in) :: A(0:NL-1)
  real(8) :: nabt(12),zx,zy,zz

  nabt( 1: 4) = nabx(1:4)
  nabt( 5: 8) = naby(1:4)
  nabt( 9:12) = nabz(1:4)

  zx = 0.d0
  zy = 0.d0
  zz = 0.d0
  call current_stencil(nabt,A,zx,zy,zz)

  zcx(ib,ik)=zx
  zcy(ib,ik)=zy
  zcz(ib,ik)=zz
end subroutine
#endif

subroutine current0
  use Global_Variables, only: NBoccmax,zu_t
  implicit none
  call current('ZE',NBoccmax,zu_t)
end subroutine

subroutine current_GS
  use Global_Variables, only: NB,zu_GS
  implicit none
  call current('GS',NB,zu_GS)
end subroutine

subroutine current_RT
  use Global_Variables, only: NBoccmax,zu_t
  implicit none
  call current('RT',NBoccmax,zu_t)
end subroutine

subroutine current_RT_MS(ixy_m)
  use Global_Variables, only: NBoccmax,zu_m
  implicit none
  integer,intent(in) :: ixy_m
  call current('RT',NBoccmax,zu_m(:,:,:,ixy_m))
end subroutine

subroutine current(mode,NBtmp,zutmp)
  use Global_Variables, only: NL,NK_s,NK_e
  use timer
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  implicit none
  character(2), intent(in) :: mode
  integer,intent(in)       :: NBtmp
  complex(8),intent(in)    :: zutmp(NL,NBtmp,NK_s:NK_e)
  real(8) :: jx,jy,jz

  NVTX_BEG('current_RT()',2)
  call timer_begin(LOG_CURRENT)
#ifdef _OPENACC
  if (mode == 'RT') then
    call current_acc_impl(zutmp,jx,jy,jz)
  else
    call impl(mode,NBtmp,zutmp,jx,jy,jz)
  end if
#else
  call impl(mode,NBtmp,zutmp,jx,jy,jz)
#endif
  call summation(jx,jy,jz)
  call timer_end(LOG_CURRENT)
  NVTX_END()

contains
  subroutine impl(mode, NBtmp, zutmp, jxs, jys, jzs)
    use Global_Variables
    implicit none
    character(2),intent(in) :: mode
    integer,intent(in)      :: NBtmp
    complex(8),intent(in)   :: zutmp(0:NL-1,NBtmp,NK_s:NK_e)
    real(8),intent(out)     :: jxs,jys,jzs

    integer :: ikb,ib,ik,i,j,ix,iy,iz,ia
    real(8) :: kr,jx,jy,jz,IaLxyz
    real(8) :: nabt(12)

    nabt( 1: 4) = nabx(1:4)
    nabt( 5: 8) = naby(1:4)
    nabt( 9:12) = nabz(1:4)

    IaLxyz = 1.0 / aLxyz

    jxs=0.d0
    jys=0.d0
    jzs=0.d0
!$omp parallel reduction(+:jxs,jys,jzs)

    !Constructing nonlocal part
!$omp do private(ik,ia,j,i,ix,iy,iz,kr) collapse(2)
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

    if (mode == 'ZE' .or. mode == 'GS') then
!$omp do private(ikb,ik,ib,jx,jy,jz)
      do ikb=1,NKB
        ik=ik_table(ikb)
        ib=ib_table(ikb)

        call init(ik,ib,zutmp(:,ib,ik),jx,jy,jz)
        call stencil(ik,ib,zutmp(:,ib,ik),nabt,jx,jy,jz)
        call pseudo_pt(ik,ib,zutmp(:,ib,ik),IaLxyz,jx,jy,jz)

        jxs=jxs+jx
        jys=jys+jy
        jzs=jzs+jz
      end do
!$omp end do
    else if (mode == 'RT') then
!$omp do private(ikb,ik,ib,jx,jy,jz)
      do ikb=1,NKB
        ik=ik_table(ikb)
        ib=ib_table(ikb)

        call init(ik,ib,zutmp(:,ib,ik),jx,jy,jz)
#ifdef ARTED_CURRENT_PREPROCESSING
        call preconditioned_stencil(ik,ib,jx,jy,jz)
#else
        call stencil(ik,ib,zutmp(:,ib,ik),nabt,jx,jy,jz)
#endif
        call pseudo_pt(ik,ib,zutmp(:,ib,ik),IaLxyz,jx,jy,jz)

        jxs=jxs+jx
        jys=jys+jy
        jzs=jzs+jz
      end do
!$omp end do
    else
      call err_finalize('ERROR: current mode')
    end if
!$omp end parallel
  end subroutine

  subroutine init(ik,ib,zutmp,jx,jy,jz)
    use Global_Variables
    implicit none
    integer,intent(in)    :: ik,ib
    complex(8),intent(in) :: zutmp(0:NL-1)
    real(8),intent(out)   :: jx,jy,jz

    integer :: i
    real(8) :: jt

    jt=0.d0
!dir$ vector aligned
    do i=0,NL-1
      jt=jt+real(zutmp(i))**2+imag(zutmp(i))**2
    end do

    jx=occ(ib,ik)*kAc(ik,1)*jt
    jy=occ(ib,ik)*kAc(ik,2)*jt
    jz=occ(ib,ik)*kAc(ik,3)*jt
  end subroutine

  subroutine stencil(ik,ib,zutmp,nabt,jx,jy,jz)
    use Global_Variables
    use opt_variables
    implicit none
    integer,intent(in)    :: ik,ib
    complex(8),intent(in) :: zutmp(0:NL-1)
    real(8),intent(in)    :: nabt(12)
    real(8),intent(inout) :: jx,jy,jz
    real(8) :: jxt,jyt,jzt

    jxt=0.d0
    jyt=0.d0
    jzt=0.d0
    call current_stencil(nabt,zutmp,jxt,jyt,jzt)

    jx=jx+jxt*occ(ib,ik)
    jy=jy+jyt*occ(ib,ik)
    jz=jz+jzt*occ(ib,ik)
  end subroutine

#ifdef ARTED_CURRENT_PREPROCESSING
  subroutine preconditioned_stencil(ik,ib,jx,jy,jz)
    use Global_Variables
    use opt_variables
    implicit none
    integer,intent(in)    :: ik,ib
    real(8),intent(inout) :: jx,jy,jz

    jx=jx+zcx(ib,ik)*occ(ib,ik)
    jy=jy+zcy(ib,ik)*occ(ib,ik)
    jz=jz+zcz(ib,ik)*occ(ib,ik)
  end subroutine
#endif

  subroutine pseudo_pt(ik,ib,zutmp,IaLxyz,jx,jy,jz)
    use Global_Variables
    use omp_lib
    use opt_variables
    use timer
    implicit none
    integer,intent(in)    :: ik,ib
    complex(8),intent(in) :: zutmp(NL)
    real(8),intent(in)    :: IaLxyz
    real(8),intent(out)   :: jx,jy,jz

    integer    :: ilma,ia,j,i,ix,iy,iz
    real(8)    :: x,y,z
    real(8)    :: jxt,jyt,jzt
    complex(8) :: uVpsi,uVpsix,uVpsiy,uVpsiz

    jxt=0d0
    jyt=0d0
    jzt=0d0

    do ilma=1,Nlma
      ia=a_tbl(ilma)
      uVpsi=0.d0; uVpsix=0.d0; uVpsiy=0.d0; uVpsiz=0.d0
      do j=1,Mps(ia)
        i=Jxyz(j,ia)

        ix=Jxx(j,ia); x=Lx(i)*Hx-ix*aLx
        iy=Jyy(j,ia); y=Ly(i)*Hy-iy*aLy
        iz=Jzz(j,ia); z=Lz(i)*Hz-iz*aLz

        uVpsi =uVpsi +uV(j,ilma)*ekr_omp(j,ia,ik)  *zutmp(i)
        uVpsix=uVpsix+uV(j,ilma)*ekr_omp(j,ia,ik)*x*zutmp(i)
        uVpsiy=uVpsiy+uV(j,ilma)*ekr_omp(j,ia,ik)*y*zutmp(i)
        uVpsiz=uVpsiz+uV(j,ilma)*ekr_omp(j,ia,ik)*z*zutmp(i)
      end do
      uVpsi =uVpsi *Hxyz*iuV(ilma)
      uVpsix=uVpsix*Hxyz
      uVpsiy=uVpsiy*Hxyz
      uVpsiz=uVpsiz*Hxyz
      jxt=jxt+occ(ib,ik)*IaLxyz*2*imag(conjg(uVpsix)*uVpsi)
      jyt=jyt+occ(ib,ik)*IaLxyz*2*imag(conjg(uVpsiy)*uVpsi)
      jzt=jzt+occ(ib,ik)*IaLxyz*2*imag(conjg(uVpsiz)*uVpsi)
    end do

    jx=jx*Hxyz*IaLxyz+jxt
    jy=jy*Hxyz*IaLxyz+jyt
    jz=jz*Hxyz*IaLxyz+jzt
  end subroutine

  subroutine summation(jx,jy,jz)
    use Global_Variables, only: jav
    use communication
    use timer
    implicit none
    real(8),intent(in) :: jx,jy,jz
    real(8) :: jav_l(3)

    call timer_begin(LOG_ALLREDUCE)
    jav_l(1)=jx
    jav_l(2)=jy
    jav_l(3)=jz
    call comm_summation(jav_l,jav,3,proc_group(2))
    call timer_end(LOG_ALLREDUCE)
  end subroutine
end subroutine
