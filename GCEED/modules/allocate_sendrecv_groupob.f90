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
MODULE allocate_sendrecv_groupob_sub

real(8),allocatable :: srmatbox1_x(:,:,:,:,:),srmatbox1_y(:,:,:,:,:),srmatbox1_z(:,:,:,:,:)
real(8),allocatable :: srmatbox2_x(:,:,:,:,:),srmatbox2_y(:,:,:,:,:),srmatbox2_z(:,:,:,:,:)
real(8),allocatable :: srmatbox3_x(:,:,:,:,:),srmatbox3_y(:,:,:,:,:),srmatbox3_z(:,:,:,:,:)
real(8),allocatable :: srmatbox4_x(:,:,:,:,:),srmatbox4_y(:,:,:,:,:),srmatbox4_z(:,:,:,:,:)

complex(8),allocatable :: scmatbox1_x(:,:,:,:,:),scmatbox1_y(:,:,:,:,:),scmatbox1_z(:,:,:,:,:)
complex(8),allocatable :: scmatbox2_x(:,:,:,:,:),scmatbox2_y(:,:,:,:,:),scmatbox2_z(:,:,:,:,:)
complex(8),allocatable :: scmatbox3_x(:,:,:,:,:),scmatbox3_y(:,:,:,:,:),scmatbox3_z(:,:,:,:,:)
complex(8),allocatable :: scmatbox4_x(:,:,:,:,:),scmatbox4_y(:,:,:,:,:),scmatbox4_z(:,:,:,:,:)

CONTAINS
!==================================================================================================
subroutine allocate_sendrecv_groupob
use scf_data
implicit none

if(iSCFRT==1.and.icalcforce==1)then

  allocate(srmatbox1_x(Nd,mg_num(2),mg_num(3),1:iobnum,1))
  allocate(srmatbox1_y(mg_num(1),Nd,mg_num(3),1:iobnum,1))
  allocate(srmatbox1_z(mg_num(1),mg_num(2),Nd,1:iobnum,1))
  allocate(srmatbox2_x(Nd,mg_num(2),mg_num(3),1:iobnum,1))
  allocate(srmatbox2_y(mg_num(1),Nd,mg_num(3),1:iobnum,1))
  allocate(srmatbox2_z(mg_num(1),mg_num(2),Nd,1:iobnum,1))
  allocate(srmatbox3_x(Nd,mg_num(2),mg_num(3),1:iobnum,1))
  allocate(srmatbox3_y(mg_num(1),Nd,mg_num(3),1:iobnum,1))
  allocate(srmatbox3_z(mg_num(1),mg_num(2),Nd,1:iobnum,1))
  allocate(srmatbox4_x(Nd,mg_num(2),mg_num(3),1:iobnum,1))
  allocate(srmatbox4_y(mg_num(1),Nd,mg_num(3),1:iobnum,1))
  allocate(srmatbox4_z(mg_num(1),mg_num(2),Nd,1:iobnum,1))

  srmatbox1_x=0.d0
  srmatbox1_y=0.d0
  srmatbox1_z=0.d0
  srmatbox2_x=0.d0
  srmatbox2_y=0.d0
  srmatbox2_z=0.d0
  srmatbox3_x=0.d0
  srmatbox3_y=0.d0
  srmatbox3_z=0.d0
  srmatbox4_x=0.d0
  srmatbox4_y=0.d0
  srmatbox4_z=0.d0
  
else if(iSCFRT==2.and.nproc_Mxin_mul/=1)then

  allocate(scmatbox1_x(Nd,mg_num(2),mg_num(3),1:iobnum,1))
  allocate(scmatbox1_y(mg_num(1),Nd,mg_num(3),1:iobnum,1))
  allocate(scmatbox1_z(mg_num(1),mg_num(2),Nd,1:iobnum,1))
  allocate(scmatbox2_x(Nd,mg_num(2),mg_num(3),1:iobnum,1))
  allocate(scmatbox2_y(mg_num(1),Nd,mg_num(3),1:iobnum,1))
  allocate(scmatbox2_z(mg_num(1),mg_num(2),Nd,1:iobnum,1))
  allocate(scmatbox3_x(Nd,mg_num(2),mg_num(3),1:iobnum,1))
  allocate(scmatbox3_y(mg_num(1),Nd,mg_num(3),1:iobnum,1))
  allocate(scmatbox3_z(mg_num(1),mg_num(2),Nd,1:iobnum,1))
  allocate(scmatbox4_x(Nd,mg_num(2),mg_num(3),1:iobnum,1))
  allocate(scmatbox4_y(mg_num(1),Nd,mg_num(3),1:iobnum,1))
  allocate(scmatbox4_z(mg_num(1),mg_num(2),Nd,1:iobnum,1))

  scmatbox1_x=0.d0
  scmatbox1_y=0.d0
  scmatbox1_z=0.d0
  scmatbox2_x=0.d0
  scmatbox2_y=0.d0
  scmatbox2_z=0.d0
  scmatbox3_x=0.d0
  scmatbox3_y=0.d0
  scmatbox3_z=0.d0
  scmatbox4_x=0.d0
  scmatbox4_y=0.d0
  scmatbox4_z=0.d0

end if

end subroutine allocate_sendrecv_groupob
!==================================================================================================

end module allocate_sendrecv_groupob_sub
