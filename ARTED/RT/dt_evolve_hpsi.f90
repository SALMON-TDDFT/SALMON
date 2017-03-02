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
#define TIMELOG_BEG(id) call timelog_thread_begin(id)
#define TIMELOG_END(id) call timelog_thread_end(id)

#ifdef ARTED_USE_NVTX
#define NVTX_BEG(name,id)  call nvtxStartRange(name,id)
#define NVTX_END()         call nvtxEndRange()
#else
#define NVTX_BEG(name,id)
#define NVTX_END()
#endif

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#ifndef ARTED_LBLK
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine dt_evolve_hpsi(flag_current)
  use Global_Variables
  use timelog
  use omp_lib
  use opt_variables
  implicit none
  integer    :: tid
  integer    :: ikb,ik,ib,i
  integer    :: iexp
  complex(8) :: zfac(4)
  logical, intent(in) :: flag_current

  zfac(1)=(-zI*dt)
  do i=2,4
    zfac(i)=zfac(i-1)*(-zI*dt)/i
  end do

  call timelog_begin(LOG_HPSI)

!$omp parallel private(tid) shared(zfac)
!$  tid=omp_get_thread_num()

!$omp do private(ik,ib,iexp)
  do ikb=1,NKB
    ik=ik_table(ikb)
    ib=ib_table(ikb)

    call init(zhtpsi(:,4,tid),zu(:,ib,ik))
    call hpsi_omp_KB_RT(ik,zhtpsi(:,4,tid),zhtpsi(:,1,tid))
    call hpsi_omp_KB_RT(ik,zhtpsi(:,1,tid),zhtpsi(:,2,tid))
    call hpsi_omp_KB_RT(ik,zhtpsi(:,2,tid),zhtpsi(:,3,tid))
    call hpsi_omp_KB_RT(ik,zhtpsi(:,3,tid),zhtpsi(:,4,tid))
    call update(zfac,zhtpsi(:,:,tid),zu(:,ib,ik))

#ifdef ARTED_CURRENT_PREPROCESSING
    if(flag_current) call current_omp_KB_ST(ib,ik,zu(:,ib,ik))
#endif
  end do
!$omp end do
!$omp end parallel

  call timelog_end(LOG_HPSI)

contains
  subroutine init(tpsi,zu)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timelog
    implicit none
    complex(8) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1)
    complex(8) :: zu(0:NLz-1,0:NLy-1,0:NLx-1)
    integer :: ix,iy,iz

    TIMELOG_BEG(LOG_HPSI_INIT)
!dir$ vector aligned
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      tpsi(iz,iy,ix)=zu(iz,iy,ix)
    end do
    end do
    end do
    TIMELOG_END(LOG_HPSI_INIT)
  end subroutine

  subroutine update(zfac,tpsi,zu)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timelog
    implicit none
    complex(8) :: zfac(4)
    complex(8) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1,4)
    complex(8) :: zu(0:NLz-1,0:NLy-1,0:NLx-1)
    integer :: ix,iy,iz

    TIMELOG_BEG(LOG_HPSI_UPDATE)
!dir$ vector aligned
    do ix=0,NLx-1
    do iy=0,NLy-1
    do iz=0,NLz-1
      zu(iz,iy,ix)=zu(iz,iy,ix)+zfac(1)*tpsi(iz,iy,ix,1) &
      &                        +zfac(2)*tpsi(iz,iy,ix,2) &
      &                        +zfac(3)*tpsi(iz,iy,ix,3) &
      &                        +zfac(4)*tpsi(iz,iy,ix,4)
    end do
    end do
    end do
    TIMELOG_END(LOG_HPSI_UPDATE)
  end subroutine
end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
! ifndef ARTED_LBLK
#else
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine dt_evolve_hpsi(flag_current)
  use Global_Variables
  use timelog
#ifdef ARTED_USE_NVTX
  use nvtx
#endif
  use omp_lib
  use opt_variables
  implicit none
  integer    :: i,iexp
  complex(8) :: zfac(4)
  integer    :: ikb_s,ikb_e
  integer    :: ikb0,ikb1,num_ikb1
  logical, intent(in) :: flag_current

  zfac(1)=(-zI*dt)
  do i=2,4
    zfac(i)=zfac(i-1)*(-zI*dt)/i
  end do

  ! NVTX_BEG('dt_evolve_hpsi()',2)
  call timelog_begin(LOG_HPSI)

!$acc data pcopy(zu) create(zhtpsi)
  do ikb0=1,NKB, blk_nkb_hpsi
    num_ikb1 = min(blk_nkb_hpsi, NKB-ikb0+1)
    ikb_s = ikb0
    ikb_e = ikb0 + num_ikb1-1

    call init_LBLK(zhtpsi(:,:,4),zu(:,:,:), ikb_s,ikb_e)
    call hpsi_acc_KB_RT_LBLK(zhtpsi(:,:,4),zhtpsi(:,:,1), ikb_s,ikb_e)
    call hpsi_acc_KB_RT_LBLK(zhtpsi(:,:,1),zhtpsi(:,:,2), ikb_s,ikb_e)
    call hpsi_acc_KB_RT_LBLK(zhtpsi(:,:,2),zhtpsi(:,:,3), ikb_s,ikb_e)
    call hpsi_acc_KB_RT_LBLK(zhtpsi(:,:,3),zhtpsi(:,:,4), ikb_s,ikb_e)
    call update_LBLK(zfac,zhtpsi(:,:,:),zu(:,:,:), ikb_s,ikb_e)

#ifdef ARTED_CURRENT_PREPROCESSING
    if(flag_current) call current_acc_KB_ST_LBLK(zu(:,:,:), ikb_s,ikb_e)
#endif
  end do
!$acc end data

  call timelog_end(LOG_HPSI)
  ! NVTX_END()

contains
  subroutine init_LBLK(tpsi,zu, ikb_s,ikb_e)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz
    use timelog
    implicit none
    integer :: ikb_s,ikb_e
    complex(8),intent(out) :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1, ikb_s:ikb_e)
    complex(8),intent(in)  :: zu(0:NLz-1,0:NLy-1,0:NLx-1, NBoccmax, NK_s:NK_e)
    integer :: ikb,ik,ib, ix,iy,iz

    TIMELOG_BEG(LOG_HPSI_INIT)
!$acc kernels pcopy(tpsi) pcopyin(zu,ib_table,ik_table)
!$acc loop gang vector(1)
    do ikb=ikb_s,ikb_e
!$acc loop collapse(3) gang vector(128)
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        ik=ik_table(ikb)
        ib=ib_table(ikb)
        tpsi(iz,iy,ix, ikb)=zu(iz,iy,ix, ib,ik)
      end do
      end do
      end do
    end do
!$acc end kernels
    TIMELOG_END(LOG_HPSI_INIT)
  end subroutine

  subroutine update_LBLK(zfac,tpsi,zu, ikb_s,ikb_e)
    use Global_Variables, only: NLx,NLy,NLz
    use opt_variables, only: PNLx,PNLy,PNLz, blk_nkb_hpsi
    use timelog
    implicit none
    integer :: ikb_s,ikb_e
    complex(8),intent(in)    :: zfac(4)
    complex(8),intent(in)    :: tpsi(0:PNLz-1,0:PNLy-1,0:PNLx-1, 0:blk_nkb_hpsi-1, 4)
    complex(8),intent(inout) :: zu(0:NLz-1,0:NLy-1,0:NLx-1, NBoccmax, NK_s:NK_e)
    integer :: ikb,ik,ib, ix,iy,iz

    TIMELOG_BEG(LOG_HPSI_UPDATE)
!$acc kernels pcopy(zu) pcopyin(tpsi,zfac,ib_table,ik_table)
!$acc loop independent gang vector(1)
    do ikb=ikb_s,ikb_e
!$acc loop collapse(3) gang vector(128)
      do ix=0,NLx-1
      do iy=0,NLy-1
      do iz=0,NLz-1
        ik=ik_table(ikb)
        ib=ib_table(ikb)
        zu(iz,iy,ix, ib,ik)=zu(iz,iy,ix, ib,ik) &
          &                +zfac(1)*tpsi(iz,iy,ix, ikb-ikb_s, 1) &
          &                +zfac(2)*tpsi(iz,iy,ix, ikb-ikb_s, 2) &
          &                +zfac(3)*tpsi(iz,iy,ix, ikb-ikb_s, 3) &
          &                +zfac(4)*tpsi(iz,iy,ix, ikb-ikb_s, 4)
      end do
      end do
      end do
    end do
!$acc end kernels
    TIMELOG_END(LOG_HPSI_UPDATE)
  end subroutine
end subroutine
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
#endif
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
