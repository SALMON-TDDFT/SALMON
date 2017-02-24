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
!This file is "hpsi.f90"
!This file contain a subroutine.
!Subroutine hpsi(q)
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
subroutine hpsi_omp_KB(ik,tpsi,ttpsi,htpsi)
  use Global_Variables, only: functional, NL
  implicit none
  integer,intent(in)     :: ik
  complex(8),intent(in)  :: tpsi(NL)
  complex(8),intent(out) :: ttpsi(NL),htpsi(NL)

  select case(functional)
    case('PZ', 'PZM','PBE','TBmBJ')
      call hpsi1(ik,tpsi,ttpsi,htpsi)
    case('TPSS','VS98')
      call err_finalize('hpsi_omp_KB: TPSS/VS98 ver. not implemented.')
  end select

contains
  subroutine hpsi1(ik,tpsi_,ttpsi_,htpsi_)
    use Global_Variables
    use opt_variables
    use timelog
    implicit none
    integer,intent(in)     :: ik
    complex(8),intent(in)  :: tpsi_(NL)
    complex(8),intent(out) :: ttpsi_(NL),htpsi_(NL)

    integer :: i,ia,j,ilma
    real(8) :: k2
    complex(8) :: uVpsi
    real(8) :: k2lap0_2
    real(8) :: nabt(12)

    call timelog_thread_begin(LOG_HPSI)

    k2=sum(kAc(ik,:)**2)
    k2lap0_2=(k2-(lapx(0)+lapy(0)+lapz(0)))*0.5d0

    call timelog_thread_begin(LOG_HPSI_STENCIL)
    nabt( 1: 4)=kAc(ik,1)*nabx(1:4)
    nabt( 5: 8)=kAc(ik,2)*naby(1:4)
    nabt( 9:12)=kAc(ik,3)*nabz(1:4)
    call hpsi1_tuned(k2lap0_2,Vloc,lapt,nabt,tpsi_,ttpsi_,htpsi_)
    call timelog_thread_end(LOG_HPSI_STENCIL)

    call timelog_thread_begin(LOG_HPSI_PSEUDO)

    !Calculating nonlocal part
    do ilma=1,Nlma
      ia=a_tbl(ilma)
      uVpsi=0.d0
      do j=1,Mps(ia)
        i=Jxyz(j,ia)
        uVpsi=uVpsi+uV(j,ilma)*ekr_omp(j,ia,ik)*tpsi_(i)
      enddo
      uVpsi=uVpsi*Hxyz*iuV(ilma)
!dir$ ivdep
      do j=1,Mps(ia)
        i=Jxyz(j,ia)
        htpsi_(i)=htpsi_(i)+conjg(ekr_omp(j,ia,ik))*uVpsi*uV(j,ilma)
      enddo
    enddo
    call timelog_thread_end(LOG_HPSI_PSEUDO)

    call timelog_thread_end(LOG_HPSI)
  end subroutine
end subroutine hpsi_omp_KB
