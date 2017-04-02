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
!This file is "Gram_Schmidt.f90"
!This file contain one subroutine.
!Subroutine Gram_Schmidt
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Gram_Schmidt
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
End Subroutine Gram_Schmidt
