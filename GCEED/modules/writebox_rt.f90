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
MODULE writebox_rt_sub

use scf_data
use new_world_sub
implicit none 
INTERFACE writebox_rt

   MODULE PROCEDURE R_writebox_rt,C_writebox_rt

END INTERFACE

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE R_writebox_rt(matbox)
!$ use omp_lib
implicit none
integer :: i1,i2,i3
real(8) :: box
real(8) :: matbox(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
real(8) :: lmatbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))
real(8) :: lmatbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))

lmatbox=0.d0
do i1=ng_sta(1),ng_end(1)
do i2=ng_sta(2),ng_end(2)
do i3=ng_sta(3),ng_end(3)
  lmatbox(i1,i2,i3)=matbox(i1,i2,i3)
end do
end do
end do

call MPI_Allreduce(lmatbox,lmatbox2, &
             lg_num(1)*lg_num(2)*lg_num(3), &
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

if(myrank.eq.0)then
  do i1=lg_sta(1),lg_end(1)
  do i2=lg_sta(2),lg_end(2)
  do i3=lg_sta(3),lg_end(3)
    box=lmatbox2(i1,i2,i3)
    write(99) box
  end do
  end do
  end do
end if

END SUBROUTINE R_writebox_rt

!=======================================================================

SUBROUTINE C_writebox_rt(matbox)
!$ use omp_lib
implicit none
integer :: i1,i2,i3
complex(8) :: box
complex(8) :: matbox(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
complex(8) :: lmatbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))
complex(8) :: lmatbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))

lmatbox=0.d0
do i1=ng_sta(1),ng_end(1)
do i2=ng_sta(2),ng_end(2)
do i3=ng_sta(3),ng_end(3)
  lmatbox(i1,i2,i3)=matbox(i1,i2,i3)
end do
end do
end do

call MPI_Allreduce(lmatbox,lmatbox2, &
             lg_num(1)*lg_num(2)*lg_num(3), &
             MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

if(myrank.eq.0)then
  do i1=lg_sta(1),lg_end(1)
  do i2=lg_sta(2),lg_end(2)
  do i3=lg_sta(3),lg_end(3)
    box=lmatbox2(i1,i2,i3)
    write(99) box
  end do
  end do
  end do
end if

END SUBROUTINE C_writebox_rt

!=======================================================================

END MODULE writebox_rt_sub

