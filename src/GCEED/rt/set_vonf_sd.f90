!
!  Copyright 2018 SALMON developers
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

subroutine set_vonf_sd
!$ use omp_lib
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation
  use scf_data
  implicit none
  integer :: i
  integer :: ix,iy,iz,iix,iiy,iiz
  integer :: max_icell
  real(8) :: rr

  if(iperiodic==0)then
    max_icell=0
  else if(iperiodic==3)then
    max_icell=2
  end if

  vonf_sd=0.d0
  eonf_sd=0.d0

  do i=1,nump
    do iiz=-max_icell,max_icell
    do iiy=-max_icell,max_icell
    do iix=-max_icell,max_icell
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        rr=sqrt((gridcoo(ix,1)-(coop(1,i)+dble(iix*lg_num(1))*Hgs(1)))**2 &
               +(gridcoo(iy,2)-(coop(2,i)+dble(iiy*lg_num(2))*Hgs(2)))**2 &
               +(gridcoo(iz,3)-(coop(3,i)+dble(iiz*lg_num(3))*Hgs(3)))**2)
        if(rr>radp_diele)then
          vonf_sd(ix,iy,iz)=vonf_sd(ix,iy,iz)  &
                        -(vecp(1,i)*(gridcoo(ix,1)-(coop(1,i)+dble(iix*lg_num(1))*Hgs(1))) &
                         +vecp(2,i)*(gridcoo(iy,2)-(coop(2,i)+dble(iiy*lg_num(2))*Hgs(2))) &
                         +vecp(3,i)*(gridcoo(iz,3)-(coop(3,i)+dble(iiz*lg_num(3))*Hgs(3))))/rr**3
          eonf_sd(1,ix,iy,iz)=eonf_sd(1,ix,iy,iz)  &
                          +(3.d0*(vecp(1,i)*(gridcoo(ix,1)-(coop(1,i)+dble(iix*lg_num(1))*Hgs(1)))/rr &
                                 +vecp(2,i)*(gridcoo(iy,2)-(coop(2,i)+dble(iiy*lg_num(2))*Hgs(2)))/rr &
                                 +vecp(3,i)*(gridcoo(iz,3)-(coop(3,i)+dble(iiz*lg_num(3))*Hgs(3)))/rr)* &
                              (gridcoo(ix,1)-(coop(1,i)+dble(iix*lg_num(1))*Hgs(1)))/rr &
                        -vecp(1,i)) /rr**3
          eonf_sd(2,ix,iy,iz)=eonf_sd(2,ix,iy,iz)  &
                          +(3.d0*(vecp(1,i)*(gridcoo(ix,1)-(coop(1,i)+dble(iix*lg_num(1))*Hgs(1)))/rr &
                                 +vecp(2,i)*(gridcoo(iy,2)-(coop(2,i)+dble(iiy*lg_num(2))*Hgs(2)))/rr &
                                 +vecp(3,i)*(gridcoo(iz,3)-(coop(3,i)+dble(iiz*lg_num(3))*Hgs(3)))/rr)* &
                              (gridcoo(iy,2)-(coop(2,i)+dble(iiy*lg_num(2))*Hgs(2)))/rr &
                        -vecp(2,i)) /rr**3
          eonf_sd(3,ix,iy,iz)=eonf_sd(3,ix,iy,iz)  &
                          +(3.d0*(vecp(1,i)*(gridcoo(ix,1)-(coop(1,i)+dble(iix*lg_num(1))*Hgs(1)))/rr &
                                 +vecp(2,i)*(gridcoo(iy,2)-(coop(2,i)+dble(iiy*lg_num(2))*Hgs(2)))/rr &
                                 +vecp(3,i)*(gridcoo(iz,3)-(coop(3,i)+dble(iiz*lg_num(3))*Hgs(3)))/rr)* &
                              (gridcoo(iz,3)-(coop(3,i)+dble(iiz*lg_num(3))*Hgs(3)))/rr &
                          -vecp(3,i)) /rr**3
        else
          vonf_sd(ix,iy,iz)=vonf_sd(ix,iy,iz)  &
                        -(vecp(1,i)*(gridcoo(ix,1)-(coop(1,i)+dble(iix*lg_num(1))*Hgs(1))) &
                         +vecp(2,i)*(gridcoo(iy,2)-(coop(2,i)+dble(iiy*lg_num(2))*Hgs(2))) &
                         +vecp(3,i)*(gridcoo(iz,3)-(coop(3,i)+dble(iiz*lg_num(3))*Hgs(3))))/radp_diele**3
          eonf_sd(1,ix,iy,iz)=eonf_sd(1,ix,iy,iz)  &
                         +vecp(1,i)/radp_diele**3
          eonf_sd(2,ix,iy,iz)=eonf_sd(2,ix,iy,iz)  &
                         +vecp(2,i)/radp_diele**3
          eonf_sd(3,ix,iy,iz)=eonf_sd(3,ix,iy,iz)  &
                         +vecp(3,i)/radp_diele**3
        end if
      end do
      end do
      end do
    end do
    end do
    end do
  end do

  return

end subroutine set_vonf_sd

