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
MODULE change_order_sub

use scf_data
use new_world_sub
implicit none 
INTERFACE change_order

   MODULE PROCEDURE R_change_order,C_change_order

END INTERFACE

CONTAINS

!======================================================================
subroutine R_change_order(tpsi)
use salmon_parallel, only: nproc_group_grid
use salmon_communication, only: comm_summation
implicit none

real(8) :: tpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,1)

integer :: is,iss,iobsta(2),iobend(2)
integer :: iob,job,iob_myob
integer :: imin,imin_myob
integer :: icheck_corrkob
real(8) :: rbox
real(8),allocatable :: matbox1(:,:,:),matbox2(:,:,:),matbox3(:,:,:),matbox4(:,:,:)

allocate(matbox1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(matbox2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(matbox3(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(matbox4(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

if(ilsda == 0)then
  iss=1
  iobsta(1)=1
  iobend(1)=itotMST
else if(ilsda == 1)then
  iss=2
  iobsta(1)=1
  iobend(1)=MST(1)
  iobsta(2)=MST(1)+1
  iobend(2)=itotMST
end if

do is=1,iss
  do iob=iobsta(is),iobend(is)-1
    call calc_myob(iob,iob_myob)
    imin=iob
    do job=iob+1,iobend(is)
      if(esp(job,1)<esp(imin,1)) imin=job
    end do
    call calc_myob(imin,imin_myob)
    if(iob/=imin)then
      rbox=esp(iob,1)
      esp(iob,1)=esp(imin,1)
      esp(imin,1)=rbox
      matbox1=0.d0
      call check_corrkob(iob,icheck_corrkob)
      if(icheck_corrkob==1) matbox1(:,:,:)=tpsi(:,:,:,iob_myob,1)
      call comm_summation(matbox1,matbox2,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_grid)
      matbox3=0.d0
      call check_corrkob(imin,icheck_corrkob)
      if(icheck_corrkob==1) matbox3(:,:,:)=tpsi(:,:,:,imin_myob,1)
      call comm_summation(matbox3,matbox4,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_grid)
      call check_corrkob(iob,icheck_corrkob)
      if(icheck_corrkob==1) tpsi(:,:,:,iob_myob,1)=matbox4(:,:,:)
      call check_corrkob(imin,icheck_corrkob)
      if(icheck_corrkob==1) tpsi(:,:,:,imin_myob,1)=matbox2(:,:,:)
    end if
  end do
end do

deallocate(matbox1,matbox2,matbox3,matbox4)

return

end subroutine R_change_order
!======================================================================
subroutine C_change_order(tpsi)
use salmon_parallel, only: nproc_group_grid
use salmon_communication, only: comm_summation
implicit none

complex(8) :: tpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,1)

integer :: is,iss,iobsta(2),iobend(2)
integer :: iob,job,iob_myob
integer :: imin,imin_myob
integer :: icheck_corrkob
real(8) :: rbox
complex(8),allocatable :: matbox1(:,:,:),matbox2(:,:,:),matbox3(:,:,:),matbox4(:,:,:)

allocate(matbox1(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(matbox2(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(matbox3(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))
allocate(matbox4(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

if(ilsda == 0)then
  iss=1
  iobsta(1)=1
  iobend(1)=itotMST
else if(ilsda == 1)then
  iss=2
  iobsta(1)=1
  iobend(1)=MST(1)
  iobsta(2)=MST(1)+1
  iobend(2)=itotMST
end if

do is=1,iss
  do iob=iobsta(is),iobend(is)-1
    call calc_myob(iob,iob_myob)
    imin=iob
    do job=iob+1,iobend(is)
      if(esp(job,1)<esp(imin,1)) imin=job
    end do
    call calc_myob(imin,imin_myob)
    if(iob/=imin)then
      rbox=esp(iob,1)
      esp(iob,1)=esp(imin,1)
      esp(imin,1)=rbox
      matbox1=0.d0
      call check_corrkob(iob,icheck_corrkob)
      if(icheck_corrkob==1) matbox1(:,:,:)=tpsi(:,:,:,iob_myob,1)
      call comm_summation(matbox1,matbox2,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_grid)
      matbox3=0.d0
      call check_corrkob(imin,icheck_corrkob)
      if(icheck_corrkob==1) matbox3(:,:,:)=tpsi(:,:,:,imin_myob,1)
      call comm_summation(matbox3,matbox4,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_grid)
      call check_corrkob(iob,icheck_corrkob)
      if(icheck_corrkob==1) tpsi(:,:,:,iob_myob,1)=matbox4(:,:,:)
      call check_corrkob(imin,icheck_corrkob)
      if(icheck_corrkob==1) tpsi(:,:,:,imin_myob,1)=matbox2(:,:,:)
    end if
  end do
end do

deallocate(matbox1,matbox2,matbox3,matbox4)

return

end subroutine C_change_order
!======================================================================

END MODULE change_order_sub
