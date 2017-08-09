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
MODULE readbox_rt_sub

use scf_data
use new_world_sub
implicit none 
INTERFACE readbox_rt

   MODULE PROCEDURE R_readbox_rt,C_readbox_rt

END INTERFACE

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE R_readbox_rt(matbox)
use salmon_parallel, only: nproc_id_global, nproc_group_global
use salmon_communication, only: comm_is_root, comm_bcast
implicit none
integer :: ix,iy,iz,i1,i2,i3
real(8) :: box
real(8) :: matbox(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
real(8) :: lmatbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))

if(comm_is_root(nproc_id_global))then
  do i1=lg_sta(1),lg_end(1)
  do i2=lg_sta(2),lg_end(2)
  do i3=lg_sta(3),lg_end(3)
    read(98) box
    lmatbox(i1,i2,i3)=box
  end do
  end do
  end do
end if

call comm_bcast(lmatbox,nproc_group_global)

do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  matbox(ix,iy,iz)=lmatbox(ix,iy,iz)
end do
end do
end do

END SUBROUTINE R_readbox_rt

!=======================================================================

SUBROUTINE C_readbox_rt(matbox)
use salmon_parallel, only: nproc_id_global, nproc_group_global
use salmon_communication, only: comm_is_root, comm_bcast
implicit none
integer :: ix,iy,iz,i1,i2,i3
complex(8) :: box
complex(8) :: matbox(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
complex(8) :: lmatbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3))

if(comm_is_root(nproc_id_global))then
  do i1=lg_sta(1),lg_end(1)
  do i2=lg_sta(2),lg_end(2)
  do i3=lg_sta(3),lg_end(3)
    read(98) box
    lmatbox(i1,i2,i3)=box
  end do
  end do
  end do
end if

call comm_bcast(lmatbox,nproc_group_global)

do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  matbox(ix,iy,iz)=lmatbox(ix,iy,iz)
end do
end do
end do

END SUBROUTINE C_readbox_rt

!=======================================================================

END MODULE readbox_rt_sub

