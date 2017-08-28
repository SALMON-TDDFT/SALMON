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
!    call hpsi_omp_KB_base(ik,zhtpsi(:,1,tid),zhtpsi(:,2,tid),zttpsi(:,tid))
    call hpsi_test3(ik,zhtpsi(:,1,tid),zhtpsi(:,2,tid),zttpsi(:,tid))
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
!    call hpsi_omp_KB_base(ik,tpsi,htpsi)
    call hpsi_test3(ik,tpsi,htpsi)
  end subroutine

!-----------------------------------------------------------------------------------------------------------------------------------
  subroutine hpsi_test3(ik,tpsi,htpsi,ttpsi)
    use hpsi_sub, only: hpsi_C
    use timer
    use Global_Variables, only: NLx,NLy,NLz,kAc,lapx,lapy,lapz,nabx,naby,nabz,Vloc,Mps,uV,iuV,Hxyz,ekr_omp,Nlma,a_tbl,Nps,NI,Jxyz
    use opt_variables, only: lapt,PNLx,PNLy,PNLz,PNL
#ifdef ARTED_USE_NVTX
    use nvtx
#endif
    implicit none
    integer,intent(in)              :: ik
    complex(8),intent(in)           ::  tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out)          :: htpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8),intent(out),optional :: ttpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    !
    integer :: ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end &
              ,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end &
              ,is_table(1),ik_table(1),i,irank_overlap(6),icomm_overlap,icomm_pseudo
    real(8) :: lap0,lapt_wrk(4,3),nabt_wrk(4,3),kAc_wrk(1,3)
    integer,allocatable :: idx(:),idy(:),idz(:),Jxyz_wrk(:,:,:)
    real(8),allocatable :: uVu_wrk(:)

    lap0 = -(lapx(0)+lapy(0)+lapz(0))*0.5d0

    lapt_wrk(1:4,1) = lapz(1:4) ! x <--> z
    lapt_wrk(1:4,2) = lapy(1:4)
    lapt_wrk(1:4,3) = lapx(1:4) ! x <--> z

    nabt_wrk(1:4,1) = nabz(1:4) ! x <--> z
    nabt_wrk(1:4,2) = naby(1:4)
    nabt_wrk(1:4,3) = nabx(1:4) ! x <--> z

    kAc_wrk(1,1) = kAc(ik,3) ! x <--> z
    kAc_wrk(1,2) = kAc(ik,2)
    kAc_wrk(1,3) = kAc(ik,1) ! x <--> z

    ix_sta = 0
    ix_end = NLz-1 ! x <--> z
    iy_sta = 0
    iy_end = NLy-1
    iz_sta = 0
    iz_end = NLx-1 ! x <--> z

    ipx_sta = 0
    ipx_end = PNLz-1 ! x <--> z
    ipy_sta = 0
    ipy_end = PNLy-1
    ipz_sta = 0
    ipz_end = PNLx-1 ! x <--> z

    is_table(1) = 1
    ik_table(1) = 1

    allocate(idx(ix_sta-4:ix_end+4),idy(iy_sta-4:iy_end+4),idz(iz_sta-4:iz_end+4))
    do i=ix_sta-4,ix_end+4
      idx(i) = mod(NLz+i,NLz) ! x <--> z
    end do
    do i=iy_sta-4,iy_end+4
      idy(i) = mod(NLy+i,NLy)
    end do
    do i=iz_sta-4,iz_end+4
      idz(i) = mod(NLx+i,NLx) ! x <--> z
    end do

    call convert_pseudo_ARTED(Jxyz_wrk,Jxyz,NI,Nps,NLy,NLz,Nlma,uVu_wrk,iuV,Hxyz)

    call hpsi_C(tpsi,htpsi,ipx_sta,ipx_end,ipy_sta,ipy_end,ipz_sta,ipz_end,1 &
           ,Vloc,ix_sta,ix_end,iy_sta,iy_end,iz_sta,iz_end,1 &
           ,idx,idy,idz,lap0,lapt_wrk,is_table,1 &
           ,NI,Nps,Nlma,a_tbl,Mps,Jxyz_wrk,uV,uVu_wrk &
           ,1,irank_overlap,icomm_overlap,icomm_pseudo &
           ,ik_table,nabt_wrk,kAc_wrk,ekr_omp(:,:,ik:ik),ttpsi)

    deallocate(idx,idy,idz,Jxyz_wrk,uVu_wrk)

    return
  contains
    subroutine convert_pseudo_ARTED(Jxyz_new,Jxyz,NI,Nps,NLy,NLz,Nlma,uVu,iuV,Hxyz)
      integer :: NI,Nps,Jxyz(Nps,NI),NLy,NLz,Nlma,iuV(Nlma)
      integer,allocatable :: Jxyz_new(:,:,:)
      real(8) :: Hxyz
      real(8),allocatable :: uVu(:)
      !
      integer    :: ia,j,i,ix,iy,iz ! i = ix*NLy*NLz + iy*NLz + iz + 1, ix=0:NLx-1, ...
      allocate(Jxyz_new(3,Nps,NI),uVu(Nlma))
      do ia=1,NI
        do j=1,Nps
          i = Jxyz(j,ia)
          iz = mod(i-1,NLz)
          iy = mod((i-1-iz)/NLz,NLy)
          ix = (i-1-iz-iy*NLy)/(NLy*NLz)
          Jxyz_new(1,j,ia) = iz ! x <--> z
          Jxyz_new(2,j,ia) = iy
          Jxyz_new(3,j,ia) = ix ! x <--> z
        end do
      end do
      do i=1,Nlma
        uVu(i) = dble(iuV(i)) * Hxyz
      end do
      return
    end subroutine convert_pseudo_ARTED
  end subroutine hpsi_test3
!-----------------------------------------------------------------------------------------------------------------------------------

end module
