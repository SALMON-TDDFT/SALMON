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
#define LOG_BEG(id) call timer_thread_begin(id)
#define LOG_END(id) call timer_thread_end(id)

#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

subroutine hpsi_acc_KB_RT_LBLK(tpsi,htpsi, ikb_s,ikb_e)
  use Global_Variables
  use opt_variables
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  implicit none
  integer :: ikb_s,ikb_e
  complex(8),intent(in)  ::  tpsi(0:PNL-1, ikb_s:ikb_e)
  complex(8),intent(out) :: htpsi(0:PNL-1, ikb_s:ikb_e)
  integer :: ikb,ik

  NVTX_BEG('hpsi_acc_KB_RT_LBLK()', 3)
  select case(functional)
    case('PZ','PZM', 'PBE','TBmBJ')
      call hpsi1_LBLK(tpsi(:,:),htpsi(:,:), ikb_s,ikb_e)
    case('TPSS','VS98')
      call err_finalize('hpsi_acc_KB_RT_LBLK: TPSS/VS98 ver. not implemented.')
  end select
  NVTX_END()

contains
  subroutine pseudo_pt_LBLK(tpsi,htpsi, ikb_s,ikb_e)
    use Global_Variables, only: Mps,uV,iuV,Hxyz,ekr_omp,Nlma,a_tbl
#ifdef ARTED_STENCIL_PADDING
    use opt_variables, only: zJxyz => zKxyz,PNL
#else
    use opt_variables, only: zJxyz,PNL
#endif
    implicit none
    integer :: ikb_s,ikb_e
    complex(8),intent(in)  ::  tpsi(0:PNL-1, ikb_s:ikb_e)
    complex(8),intent(out) :: htpsi(0:PNL-1, ikb_s:ikb_e)

    integer    :: ilma,ia,j,i
    integer    :: ikb,ik
    complex(8) :: uVpsi0
    complex(8) :: uVpsi(Nlma, ikb_s:ikb_e)
    complex(8) :: tpsi0
    integer    :: my_nlma, n, vi

    !Calculating nonlocal part

!$acc kernels pcopy(tpsi) create(uVpsi) pcopyin(a_tbl,ekr_omp,ik_table,mps,uv,iuv) &
#ifdef ARTED_STENCIL_PADDING
!$acc& pcopyin(zKxyz) &
#else
!$acc& pcopyin(zJxyz) &
#endif
!$acc& pcopyin(t4ppt_ilma,t4ppt_j,t4ppt_nlma,t4ppt_vi2i)
!$acc loop gang vector(1)
    do ikb = ikb_s, ikb_e
!$acc loop gang vector(1)
      do ilma=1,Nlma
        ik=ik_table(ikb)
        ia=a_tbl(ilma)
        uVpsi0=cmplx(0.d0, 0.d0, kind=8)
!$acc loop vector(128) reduction(+:uVpsi0)
        do j=1,Mps(ia)
          i=zJxyz(j,ia)
          uVpsi0=uVpsi0 + uV(j,ilma)*ekr_omp(j,ia,ik)*tpsi(i,ikb)
        enddo
        uVpsi(ilma,ikb)=uVpsi0*Hxyz*iuV(ilma)
      enddo
    enddo

!$acc loop gang vector(1)
    do ikb = ikb_s, ikb_e
!$acc loop independent gang vector(128)
      do vi = 0, t4ppt_max_vi-1
        ik=ik_table(ikb)
        my_nlma = t4ppt_nlma(vi)
        if (my_nlma < 1) cycle
        tpsi0=cmplx(0.d0, 0.d0, kind=8)
!$acc loop seq
        do n = 1, my_nlma
          ilma = t4ppt_ilma(vi,n)
          j    = t4ppt_j   (vi,n)
          ia   = a_tbl(ilma)
          tpsi0= tpsi0+conjg(ekr_omp(j,ia,ik))*uVpsi(ilma,ikb)*uV(j,ilma)
        enddo
        i = t4ppt_vi2i(vi)
        htpsi(i,ikb)=htpsi(i,ikb)+tpsi0
      enddo
    enddo
!$acc end kernels

  end subroutine

  subroutine hpsi1_LBLK(tpsi,htpsi, ikb_s,ikb_e)
    use Global_Variables
    use opt_variables
    use timer
    implicit none
    integer :: ikb_s,ikb_e
    complex(8),intent(in)  ::  tpsi(0:PNL-1, ikb_s:ikb_e)
    complex(8),intent(out) :: htpsi(0:PNL-1, ikb_s:ikb_e)
    integer :: ikb,ik

    real(8) :: k2
    real(8) :: k2lap0_2(ikb_s:ikb_e)
    real(8) :: nabt(12, ikb_s:ikb_e)
!$acc data pcopy(tpsi) create(k2lap0_2,nabt)

    NVTX_BEG('hpsi1_LBLK(): hpsi1_RT_stencil', 4)
    LOG_BEG(LOG_HPSI_STENCIL)
!$acc kernels pcopy(k2lap0_2,nabt) pcopyin(ik_table,kac,lapx,lapy,lapz,nabx,naby,nabz)
!$acc loop gang vector
    do ikb = ikb_s, ikb_e
      ik=ik_table(ikb)
      k2=sum(kAc(ik,:)**2)
      k2lap0_2(ikb)=(k2-(lapx(0)+lapy(0)+lapz(0)))*0.5d0
      nabt( 1: 4,ikb)=kAc(ik,1)*nabx(1:4)
      nabt( 5: 8,ikb)=kAc(ik,2)*naby(1:4)
      nabt( 9:12,ikb)=kAc(ik,3)*nabz(1:4)
    enddo
!$acc end kernels
    call hpsi1_RT_stencil_LBLK(k2lap0_2(:),Vloc,lapt,nabt(:,:),tpsi(:,:),htpsi(:,:), ikb_s,ikb_e)
    LOG_END(LOG_HPSI_STENCIL)
    NVTX_END()

    NVTX_BEG('hpsi1_LBLK(): pseudo_pt', 5)
    LOG_BEG(LOG_HPSI_PSEUDO)
    call pseudo_pt_LBLK(tpsi(:,:),htpsi(:,:), ikb_s,ikb_e)
    LOG_END(LOG_HPSI_PSEUDO)
    NVTX_END()

!$acc end data
  end subroutine
end subroutine hpsi_acc_KB_RT_LBLK
