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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
#define LOG_BEG(id) call timer_thread_begin(id)
#define LOG_END(id) call timer_thread_end(id)

#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

module hpsi
  use timer
  implicit none

contains
  subroutine hpsi_omp_KB_GS(ik,tpsi,ttpsi,htpsi)
    use Global_Variables, only: NL,NLz,NLy,NLx
    use opt_variables, only: zhtpsi,zttpsi,PNLx,PNLy,PNLz
    use omp_lib
    implicit none
    integer,intent(in)     :: ik
    complex(8),intent(in)  :: tpsi(NL)
    complex(8),intent(out) :: ttpsi(NL),htpsi(NL)
    integer :: tid

    LOG_BEG(LOG_HPSI)

    tid = omp_get_thread_num()
    call init(tpsi,zhtpsi(:,1,tid))
    call hpsi_omp_KB_base(ik,zhtpsi(:,1,tid),zhtpsi(:,2,tid),zttpsi(:,tid))
    call copyout(zhtpsi(:,2,tid),zttpsi(:,tid),htpsi,ttpsi)

   LOG_END(LOG_HPSI)

  contains
      subroutine init(zu,tpsi)
      implicit none
      complex(8),intent(in)  :: zu(0:NLz-1,0:NLy-1,0:NLx-1)
      complex(8),intent(out) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      integer :: ix,iy,iz

!dir$ vector aligned
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        tpsi(iz,iy,ix)=zu(iz,iy,ix)
      end do
      end do
      end do
    end subroutine

    subroutine copyout(zhtpsi,zttpsi,htpsi,ttpsi)
      implicit none
      complex(8), intent(in)  :: zhtpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      complex(8), intent(in)  :: zttpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      complex(8), intent(out) :: htpsi(0:NLz-1,0:NLy-1,0:NLx-1)
      complex(8), intent(out) :: ttpsi(0:NLz-1,0:NLy-1,0:NLx-1)
      integer :: ix,iy,iz

!dir$ vector aligned
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        htpsi(iz,iy,ix) = zhtpsi(iz,iy,ix)
        ttpsi(iz,iy,ix) = zttpsi(iz,iy,ix)
      end do
      end do
      end do
    end subroutine
  end subroutine

  subroutine hpsi_omp_KB_RT(ik,tpsi,htpsi)
    use opt_variables, only: PNLx,PNLy,PNLz
    implicit none
    integer,intent(in)     :: ik
    complex(8),intent(in)  ::  tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out) :: htpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    call hpsi_omp_KB_base(ik,tpsi,htpsi)
  end subroutine

  subroutine hpsi_omp_KB_base(ik,tpsi,htpsi,ttpsi)
    use timer
    use Global_Variables, only: NLx,NLy,NLz,kAc,lapx,lapy,lapz,nabx,naby,nabz,Vloc,Mps,iuV,Hxyz,Nlma,zproj, & 
    & flag_set_ini_Ac_alocal, Ac2_al,Ac1x_al,Ac1y_al,Ac1z_al,nabt_al
    use opt_variables, only: lapt,PNLx,PNLy,PNLz,PNL,spseudo,dpseudo
#ifdef ARTED_USE_NVTX
    use nvtx
#endif
    use salmon_parallel, only: get_thread_id
    use salmon_global, only: alocal_laser
    use Ac_alocal_laser
    implicit none
    integer,intent(in)              :: ik
    complex(8),intent(in)           ::  tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out)          :: htpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out),optional :: ttpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    real(8) :: k2,k2lap0_2
    real(8) :: nabt(12)
    integer :: tid

    NVTX_BEG('hpsi1()',3)

    k2=sum(kAc(ik,:)**2)
    k2lap0_2=(k2-(lapx(0)+lapy(0)+lapz(0)))*0.5d0
    nabt( 1: 4)=kAc(ik,1)*nabx(1:4)
    nabt( 5: 8)=kAc(ik,2)*naby(1:4)
    nabt( 9:12)=kAc(ik,3)*nabz(1:4)

    LOG_BEG(LOG_HPSI_STENCIL)
      call hpsi1_RT_stencil(k2lap0_2,Vloc,lapt,nabt,tpsi,htpsi)
      if(alocal_laser=='y' .and. flag_set_ini_Ac_alocal)then
         call hpsi1_RT_stencil_add_Ac_alocal(Ac2_al(:,ik),Ac1x_al,Ac1y_al,Ac1z_al,nabt_al,tpsi,htpsi)
      endif
      if (present(ttpsi)) then
        call subtraction(Vloc,tpsi,htpsi,ttpsi)
      end if
    LOG_END(LOG_HPSI_STENCIL)

    LOG_BEG(LOG_HPSI_PSEUDO)
      tid = get_thread_id()
      call pseudo_pt(ik,zproj(:,:,ik),tpsi,htpsi,spseudo(:,tid),dpseudo(:,tid))
    LOG_END(LOG_HPSI_PSEUDO)

    NVTX_END()

  contains
    subroutine subtraction(Vloc,tpsi,htpsi,ttpsi)
      implicit none
      real(8),    intent(in)  :: Vloc(0:NLz-1,0:NLy-1,0:NLx-1)
      complex(8), intent(in)  ::  tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      complex(8), intent(in)  :: htpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      complex(8), intent(out) :: ttpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
      integer :: ix,iy,iz

!dir$ vector aligned
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        ttpsi(iz,iy,ix) = htpsi(iz,iy,ix) - Vloc(iz,iy,ix)*tpsi(iz,iy,ix)
      end do
      end do
      end do
    end subroutine

    !Calculating nonlocal part
    subroutine pseudo_pt(ik,zproj,tpsi,htpsi,spseudo,dpseudo)
      use opt_variables, only: NPI,pseudo_start_idx,idx_proj,nprojector,idx_lma
      use global_variables, only: NI,Nps
      implicit none
      integer,    intent(in)  :: ik
      complex(8), intent(in)  :: zproj(Nps,Nlma)
      complex(8), intent(in)  :: tpsi(0:PNL-1)
      complex(8), intent(out) :: htpsi(0:PNL-1)
      complex(8)              :: spseudo(NPI),dpseudo(NPI) ! working mem.

      integer    :: ia,i,j,ip,ioffset
      complex(8) :: uVpsi

      dpseudo = cmplx(0.d0)

      ! gather (load) pseudo potential point
      do i=1,NPI
        spseudo(i) = tpsi(idx_proj(i))
      end do

      do ia=1,NI
      do ip=1,nprojector(ia)
        i = idx_lma(ia) + ip

        ! summarize vector
        uVpsi   = 0.d0
        ioffset = pseudo_start_idx(ia)
        do j=1,Mps(ia)
          uVpsi = uVpsi + conjg(zproj(j,i)) * spseudo(ioffset+j)
        end do
        uVpsi = uVpsi * Hxyz * iuV(i)

        ! apply vector
        ioffset = pseudo_start_idx(ia)
        do j=1,Mps(ia)
          dpseudo(ioffset+j) = dpseudo(ioffset+j) + zproj(j,i) * uVpsi
        end do
      end do
      end do

      ! scatter (store) pseudo potential point
      do i=1,NPI
        htpsi(idx_proj(i)) = htpsi(idx_proj(i)) + dpseudo(i)
      end do
    end subroutine
  end subroutine
end module

