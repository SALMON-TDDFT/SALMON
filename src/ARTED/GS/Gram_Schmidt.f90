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
!This file is "Gram_Schmidt.f90"
!This file contain one subroutine.
!Subroutine Gram_Schmidt
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
Subroutine Gram_Schmidt
  use Global_Variables

  select case (omp_loop)
  case('k')
     call Gram_Schmidt_ompk
  case('b')
     call Gram_Schmidt_ompb
  end select

End Subroutine Gram_Schmidt

Subroutine Gram_Schmidt_ompk
  use Global_Variables
  use timer
  implicit none
  integer :: ik,ib,ibt
  real(8) :: s
  complex(8) :: zov

  call timer_begin(LOG_GRAM_SCHMIDT)
!$omp parallel do private(ib,ibt,zov,s)
  do ik=NK_s,NK_e
  do ib=1,NB
    do ibt=1,ib-1
      zov=sum(conjg(zu_GS(:,ibt,ik))*zu_GS(:,ib,ik))*Hxyz
      zu_GS(:,ib,ik)=zu_GS(:,ib,ik)-zu_GS(:,ibt,ik)*zov
    enddo
    s=sum(abs(zu_GS(:,ib,ik))**2)*Hxyz
    zu_GS(:,ib,ik)=zu_GS(:,ib,ik)/sqrt(s)
  enddo
  enddo
  call timer_end(LOG_GRAM_SCHMIDT)

  return
End Subroutine Gram_Schmidt_ompk

Subroutine Gram_Schmidt_ompb
  !This subroutine follows the algorithm developed 
  !in RSDFT for efficient parallelization
  !(but not perfectly the same to avoid complication.
  ! If further tuning is necessary, you can revise refering the algorithm.
  use Global_Variables
  use timer
  implicit none
  integer :: i,ik,ib,ibt
  real(8) :: s1,s2
  complex(8) :: zov,zov1a,zov1b,zov2a,zov2b
  complex(8),allocatable :: zu_GSold(:,:,:),zov_omp(:,:,:)
  integer :: thr_id,omp_get_thread_num
  thr_id=0

  call timer_begin(LOG_GRAM_SCHMIDT)
  allocate(zu_GSold(NL,NB,NK_s:NK_e))
  allocate(zov_omp(2,NB,0:NUMBER_THREADS-1))


  do ik=NK_s,NK_e

  zu_GSold(:,:,ik) = zu_GS(:,:,ik)

  do ib=1,NB-1,2
     zov = sum(conjg(zu_GS(:,ib,ik))*zu_GSold(:,ib+1,ik))*Hxyz
     zu_GS(:,ib+1,ik) = zu_GS(:,ib+1,ik) - zu_GS(:,ib,ik)*zov

     s1=0d0
     s2=0d0
!$omp parallel do private(i) reduction(+:s1,s2)
     do i=1,NL
        s1 = s1 + abs(zu_GS(i,ib,  ik))**2
        s2 = s2 + abs(zu_GS(i,ib+1,ik))**2
     enddo
!$omp end parallel do
     s1 = 1d0/sqrt(s1*Hxyz)
     s2 = 1d0/sqrt(s2*Hxyz)

     zu_GS(:,ib,  ik) = zu_GS(:,ib,  ik)*s1
     zu_GS(:,ib+1,ik) = zu_GS(:,ib+1,ik)*s2

     zov_omp(:,:,:)=0d0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ibt,i) collapse(2)
     do ibt=ib+2,NB-1,2
        do i=1,NL
            zov_omp(1,ibt,  thr_id) = zov_omp(1,ibt,  thr_id)+ conjg(zu_GS(i,ib,  ik)) * zu_GS(i,ibt,  ik)
            zov_omp(2,ibt,  thr_id) = zov_omp(2,ibt,  thr_id)+ conjg(zu_GS(i,ib+1,ik)) * zu_GS(i,ibt,  ik)
            zov_omp(1,ibt+1,thr_id) = zov_omp(1,ibt+1,thr_id)+ conjg(zu_GS(i,ib,  ik)) * zu_GS(i,ibt+1,ik)
            zov_omp(2,ibt+1,thr_id) = zov_omp(2,ibt+1,thr_id)+ conjg(zu_GS(i,ib+1,ik)) * zu_GS(i,ibt+1,ik)
        enddo
     enddo
!$omp end do

     do ibt=ib+2,NB-1,2
        zov1a = sum(zov_omp(1,ibt,  :))*Hxyz
        zov1b = sum(zov_omp(2,ibt,  :))*Hxyz
        zov2a = sum(zov_omp(1,ibt+1,:))*Hxyz
        zov2b = sum(zov_omp(2,ibt+1,:))*Hxyz

!$omp do private(i)
        do i=1,NL
           zu_GS(i,ibt,  ik)= zu_GS(i,ibt,  ik)- zu_GS(i,ib,ik)*zov1a-zu_GS(i,ib+1,ik)*zov1b
           zu_GS(i,ibt+1,ik)= zu_GS(i,ibt+1,ik)- zu_GS(i,ib,ik)*zov2a-zu_GS(i,ib+1,ik)*zov2b
        enddo
!$omp end do

     enddo
!$omp end parallel

  enddo

  enddo !ik


  deallocate(zov_omp)
  deallocate(zu_GSold)
  call timer_end(LOG_GRAM_SCHMIDT)

  return
End Subroutine Gram_Schmidt_ompb


