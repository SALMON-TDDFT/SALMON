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
subroutine calcEstatic
!$use omp_lib
use scf_data
use new_world_sub
use sendrecvh_sub
implicit none
integer :: ist,ix,iy,iz
complex(8),parameter :: zi=(0.d0,1.d0)
real(8) :: Vh_wk(ng_sta(1)-Ndh:ng_end(1)+Ndh,   &
                 ng_sta(2)-Ndh:ng_end(2)+Ndh,   &
                 ng_sta(3)-Ndh:ng_end(3)+Ndh)
complex(8) :: Ex_static2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
complex(8) :: Ey_static2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
complex(8) :: Ez_static2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))

iwk_size=11
call make_iwksta_iwkend


!$OMP parallel do
do iz=ng_sta(3)-Ndh,ng_end(3)+Ndh
do iy=ng_sta(2)-Ndh,ng_end(2)+Ndh
do ix=ng_sta(1)-Ndh,ng_end(1)+Ndh
  Vh_wk(ix,iy,iz)=0.d0
end do
end do
end do

if(mod(itt,2)==1)then
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
!    Vh_wk(ix,iy,iz) = Vh_stock2(ix,iy,iz)-Vh0(ix,iy,iz)
    Vh_wk(ix,iy,iz) = Vh_stock2(ix,iy,iz)
  end do
  end do
  end do
else
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
!    Vh_wk(ix,iy,iz) = Vh_stock1(ix,iy,iz)-Vh0(ix,iy,iz)
    Vh_wk(ix,iy,iz) = Vh_stock1(ix,iy,iz)
  end do
  end do
  end do
end if

call sendrecvh(Vh_wk)

if(ng_sta(1)==lg_sta(1))then
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
    do ix=1,Ndh
      Vh_wk(ng_sta(1)-ix,iy,iz) = Vh_wk(ng_sta(1),iy,iz)
    end do
  end do
  end do
end if

if(ng_end(1)==lg_end(1))then
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
    do ix=1,Ndh
      Vh_wk(ng_end(1)+ix,iy,iz) = Vh_wk(ng_end(1),iy,iz)
    end do
  end do
  end do
end if

if(ng_sta(2)==lg_sta(2))then
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do ix=ng_sta(1),ng_end(1)
    do iy=1,Ndh
      Vh_wk(ix,ng_sta(2)-iy,iz) = Vh_wk(ix,ng_sta(2),iz)
    end do
  end do
  end do
end if

if(ng_end(2)==lg_end(2))then
!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do ix=ng_sta(1),ng_end(1)
    do iy=1,Ndh
      Vh_wk(ix,ng_end(2)+iy,iz) = Vh_wk(ix,ng_end(2),iz)
    end do
  end do
  end do
end if

if(ng_sta(3)==lg_sta(3))then
!$OMP parallel do
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    do iz=1,Ndh
      Vh_wk(ix,iy,ng_sta(3)-iz) = Vh_wk(ix,iy,ng_sta(3))
    end do
  end do
  end do
end if

if(ng_end(3)==lg_end(3))then
!$OMP parallel do
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    do iz=1,Ndh
      Vh_wk(ix,iy,ng_end(3)+iz) = Vh_wk(ix,iy,ng_end(3))
    end do
  end do
  end do
end if

!$OMP parallel do
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  Ex_static2(ix,iy,iz)=0.d0
  Ey_static2(ix,iy,iz)=0.d0
  Ez_static2(ix,iy,iz)=0.d0
end do
end do
end do

!$OMP parallel do
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  do ist=1,Ndh
    Ex_static2(ix,iy,iz)=Ex_static2(ix,iy,iz)               &
       -(bNmat(ist,Ndh)*(Vh_wk(ix+ist,iy,iz)-Vh_wk(ix-ist,iy,iz)))/Hgs(1)
    Ey_static2(ix,iy,iz)=Ey_static2(ix,iy,iz)               &
       -(bNmat(ist,Ndh)*(Vh_wk(ix,iy+ist,iz)-Vh_wk(ix,iy-ist,iz)))/Hgs(2)
    Ez_static2(ix,iy,iz)=Ez_static2(ix,iy,iz)               &
       -(bNmat(ist,Ndh)*(Vh_wk(ix,iy,iz+ist)-Vh_wk(ix,iy,iz-ist)))/Hgs(3)
  end do
end do
end do
end do

call MPI_Allreduce(Ex_static2,Ex_static,mg_num(1)*mg_num(2)*mg_num(3), &
             MPI_DOUBLE_COMPLEX,MPI_SUM,newworld_comm_grid,ierr)
call MPI_Allreduce(Ey_static2,Ey_static,mg_num(1)*mg_num(2)*mg_num(3), &
             MPI_DOUBLE_COMPLEX,MPI_SUM,newworld_comm_grid,ierr)
call MPI_Allreduce(Ez_static2,Ez_static,mg_num(1)*mg_num(2)*mg_num(3), &
             MPI_DOUBLE_COMPLEX,MPI_SUM,newworld_comm_grid,ierr)

end subroutine calcEstatic

