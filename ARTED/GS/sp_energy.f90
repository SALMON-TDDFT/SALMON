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
!This file is "sp_energy.f90"
!This file contain one subroutine.
!SUBROUTINE sp_energy
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine sp_energy_omp
  use Global_Variables
  use communication
  use timelog
  implicit none
  integer :: ik,ib
  real(8) :: esp_l(NB,NK)
! sato
  integer :: ia,j,i,ix,iy,iz
  real(8) :: kr
! omp
  integer :: thr_id,omp_get_thread_num

  call timelog_begin(LOG_SP_ENERGY)
  esp_l=0.d0
  thr_id=0
!$omp parallel private(thr_id)
!$  thr_id=omp_get_thread_num()
!$omp do private(ik,ia,j,i,ix,iy,iz,kr,ib)
  do ik=NK_s,NK_e
!Constructing nonlocal part ! sato
  do ia=1,NI
    do j=1,Mps(ia)
      i=Jxyz(j,ia); ix=Jxx(j,ia); iy=Jyy(j,ia); iz=Jzz(j,ia)
      kr=kAc(ik,1)*(Lx(i)*Hx-ix*aLx)+kAc(ik,2)*(Ly(i)*Hy-iy*aLy)+kAc(ik,3)*(Lz(i)*Hz-iz*aLz)
      ekr_omp(j,ia,ik)=exp(zI*kr)
    enddo
  enddo
    do ib=1,NB
      tpsi_omp(1:NL,thr_id)=zu_GS(1:NL,ib,ik)
      call hpsi_omp_KB_GS(ik,tpsi_omp(:,thr_id),ttpsi_omp(:,thr_id),htpsi_omp(:,thr_id))
      esp_l(ib,ik)=sum(conjg(zu_GS(:,ib,ik))*htpsi_omp(:,thr_id))*Hxyz
    enddo
  enddo

!$omp end parallel

  call timelog_begin(LOG_ALLREDUCE)
  call comm_summation(esp_l,esp,NB*NK,proc_group(2))
  call timelog_end(LOG_ALLREDUCE)

  call timelog_end(LOG_SP_ENERGY)

  return
End Subroutine sp_energy_omp
