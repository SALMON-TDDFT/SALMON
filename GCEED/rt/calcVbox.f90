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
!=======================================================================
!=======================================================================

SUBROUTINE calcVbox
!$ use omp_lib
use scf_data

implicit none
integer :: ix,iy,iz,jj
integer :: ix_sta_Vbox(3),ix_end_Vbox(3)

elp3(511)=MPI_Wtime()

if(icalcforce==1.or.iflag_md==1)then
  do jj=1,3
    if(lg_sta(jj)==mg_sta(jj))then
      ix_sta_Vbox(jj)=mg_sta(jj)
    else
      ix_sta_Vbox(jj)=mg_sta(jj)-Nd
    end if
    if(lg_end(jj)==mg_end(jj))then
      ix_end_Vbox(jj)=mg_end(jj)
    else
      ix_end_Vbox(jj)=mg_end(jj)+Nd
    end if
  end do
else
  ix_sta_Vbox(1:3)=mg_sta(1:3)
  ix_end_Vbox(1:3)=mg_end(1:3)
end if

select case(ikind_eext)
  case(1,2,6)
    if(dt*itt <= tau)then
!$OMP parallel do private(ix,iy,iz)
      do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
      do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
      do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
        Vbox(ix,iy,iz)=Veff(ix,iy,iz)*sin(romega*dble(itt)*dt)*sin(Pi*dble(itt)*dt/pulse_T)**2 
      end do
      end do
      end do
    else
!$OMP parallel do private(ix,iy,iz)
      do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
      do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
      do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
        Vbox(ix,iy,iz)=0.d0
      end do
      end do
      end do
    end if
  case(4)
    if(dt*itt <= tau2(1)) then
!$OMP parallel do private(ix,iy,iz)
      do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
      do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
      do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
        Vbox(ix,iy,iz)=Veff2(ix,iy,iz,1)      &
                    *sin(Pi*dt*itt/pulse_T2(1))**2      &
                    *cos(romega2(1)*dt*itt)      &
                +Veff2(ix,iy,iz,2)      &
                    *sin(Pi*dt*itt/pulse_T2(2))**2      &
                    *sin(romega2(2)*dt*itt)      
      end do
      end do
      end do
    else
!$OMP parallel do private(ix,iy,iz)
      do iz=ix_sta_Vbox(3),ix_end_Vbox(3)
      do iy=ix_sta_Vbox(2),ix_end_Vbox(2)
      do ix=ix_sta_Vbox(1),ix_end_Vbox(1)
        Vbox(ix,iy,iz)=0.d0
      end do
      end do
      end do
    endif

end select

return

END SUBROUTINE calcVbox
