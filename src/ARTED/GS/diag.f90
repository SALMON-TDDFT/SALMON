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
!This file is "diag.f"
!This file contain two subroutines.
!SUBROUTINE diag
!SUBROUTINE hermit_jacobi(n,ar,ai,vr,vi,er)

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine diag_omp
  use Global_Variables
  use timer
  use omp_lib
  use hpsi, only: hpsi_omp_KB_GS
  use projector
  implicit none
  integer,parameter :: matz=1
  integer           :: ik,ib1,ib2
  integer           :: thr_id

!LAPACK
  integer                :: lwork,info
  complex(8),allocatable :: za(:,:)
  complex(8),allocatable :: work_lp(:)
  real(8),allocatable    :: rwork(:),w(:)
#ifdef __FUJITSU
  complex(8),allocatable :: zutmp(:,:,:)
  integer                :: ithr
#else
  complex(8),allocatable :: zutmp(:,:)
#endif

  lwork=6*NB
  thr_id=0

  call timer_begin(LOG_DIAG)
  call update_projector(kac)

! FIXME: For Fujitsu's parallelized BLAS routines, the application crashes when
!        invoking a routine under OpenMP parallelized loop.
!        However, for many-core processors, outermost loop must be parallelized
!        due to OpenMP overheads.
#ifdef __FUJITSU
  allocate(za(NB,NB),zutmp(NL,NB,0:NUMBER_THREADS-1))
  allocate(work_lp(lwork),rwork(3*NB-2),w(NB))
  do ik=NK_s,NK_e
!$omp parallel private(thr_id)
!$ thr_id = omp_get_thread_num()
!$omp do private(ib1,ib2)
    do ib1=1,NB
      tpsi_omp(1:NL,thr_id)=zu_GS(1:NL,ib1,ik)
      call hpsi_omp_KB_GS(ik,tpsi_omp(:,thr_id),ttpsi_omp(:,thr_id),htpsi_omp(:,thr_id))
      do ib2=ib1+1,NB
        za(ib2,ib1)=sum(conjg(zu_GS(:,ib2,ik))*htpsi_omp(:,thr_id))*Hxyz
        za(ib1,ib2)=conjg(za(ib2,ib1))
      end do
      za(ib1,ib1)=real(sum(conjg(zu_GS(:,ib1,ik))*htpsi_omp(:,thr_id))*Hxyz)
    end do
!$omp end do
!$omp end parallel
    call zheev('V', 'U', NB, za, NB, w, work_lp, lwork, rwork, info)

    zutmp=0.d0
!$omp parallel private(thr_id)
!$ thr_id = omp_get_thread_num()
!$omp do private(ib1,ib2) collapse(2)
    do ib1=1,NB
    do ib2=1,NB
      zutmp(:,ib1,thr_id)=zutmp(:,ib1,thr_id)+zu_GS(:,ib2,ik)*za(ib2,ib1)
    end do
    end do
!$omp end parallel

    Zu_GS(:,:,ik)=0.d0
    do ithr=0,NUMBER_THREADS-1
      zu_GS(:,:,ik)=Zu_GS(:,:,ik)+zutmp(:,:,ithr)
    end do
    esp(:,ik)=w(:)
  enddo
  deallocate(za,zutmp,work_lp,rwork)
#else


!  if( NK_e-NK_s+1 .lt. NB ) then !change as you like

  select case (omp_loop)
  case('k')

  !!---Openmp is for k-points---

!$omp parallel private(thr_id,za,zutmp,work_lp,rwork,w)
!$ thr_id = omp_get_thread_num()
  allocate(za(NB,NB),zutmp(NL,NB))
  allocate(work_lp(lwork),rwork(3*NB-2),w(NB))
!$omp do private(ik,ib1,ib2)
  do ik=NK_s,NK_e
    do ib1=1,NB
      tpsi_omp(1:NL,thr_id)=zu_GS(1:NL,ib1,ik)
      call hpsi_omp_KB_GS(ik,tpsi_omp(:,thr_id),ttpsi_omp(:,thr_id),htpsi_omp(:,thr_id))
      do ib2=ib1+1,NB
        za(ib2,ib1)=sum(conjg(zu_GS(:,ib2,ik))*htpsi_omp(:,thr_id))*Hxyz
        za(ib1,ib2)=conjg(za(ib2,ib1))
      end do
      za(ib1,ib1)=real(sum(conjg(zu_GS(:,ib1,ik))*htpsi_omp(:,thr_id))*Hxyz)
    end do
    call zheev('V', 'U', NB, za, NB, w, work_lp, lwork, rwork, info)

    zutmp=0.d0
    do ib1=1,NB
    do ib2=1,NB
      zutmp(:,ib1)=zutmp(:,ib1)+zu_GS(:,ib2,ik)*za(ib2,ib1)
    end do
    end do
    zu_GS(:,:,ik)=zutmp(:,:)
    esp(:,ik)=w(:)
  enddo
!$omp end do
  deallocate(za,zutmp,work_lp,rwork)
!$omp end parallel

  case('b')

  !!---Openmp is for band orbital---

  allocate(za(NB,NB),zutmp(NL,NB))
  allocate(work_lp(lwork),rwork(3*NB-2),w(NB))
  do ik=NK_s,NK_e
!$omp parallel private(thr_id)
!$ thr_id = omp_get_thread_num()
!$omp do private(ib1,ib2)
    do ib1=1,NB
      tpsi_omp(1:NL,thr_id)=zu_GS(1:NL,ib1,ik)
      call hpsi_omp_KB_GS(ik,tpsi_omp(:,thr_id),ttpsi_omp(:,thr_id),htpsi_omp(:,thr_id))
      do ib2=ib1+1,NB
        za(ib2,ib1)=sum(conjg(zu_GS(:,ib2,ik))*htpsi_omp(:,thr_id))*Hxyz
        za(ib1,ib2)=conjg(za(ib2,ib1))
      end do
      za(ib1,ib1)=real(sum(conjg(zu_GS(:,ib1,ik))*htpsi_omp(:,thr_id))*Hxyz)
    end do
!$omp end do
!$omp end parallel
    call zheev('V', 'U', NB, za, NB, w, work_lp, lwork, rwork, info)

    zutmp=0.d0
    do ib1=1,NB
    do ib2=1,NB
      zutmp(:,ib1)=zutmp(:,ib1)+zu_GS(:,ib2,ik)*za(ib2,ib1)
    end do
    end do
    zu_GS(:,:,ik)=zutmp(:,:)
    esp(:,ik)=w(:)
  enddo
  deallocate(za,zutmp,work_lp,rwork)

  end select

#endif
  call timer_end(LOG_DIAG)

  return
End Subroutine diag_omp
