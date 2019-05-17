!
!  Copyright 2017-2019 SALMON developers
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
subroutine allocate_sendrecv
  use scf_data
  implicit none
  
  allocate(srmatbox1_x_3d(Nd,mg_num(2),mg_num(3)))
  allocate(srmatbox1_y_3d(mg_num(1),Nd,mg_num(3)))
  allocate(srmatbox1_z_3d(mg_num(1),mg_num(2),Nd))
  allocate(srmatbox2_x_3d(Nd,mg_num(2),mg_num(3)))
  allocate(srmatbox2_y_3d(mg_num(1),Nd,mg_num(3)))
  allocate(srmatbox2_z_3d(mg_num(1),mg_num(2),Nd))
  allocate(srmatbox3_x_3d(Nd,mg_num(2),mg_num(3)))
  allocate(srmatbox3_y_3d(mg_num(1),Nd,mg_num(3)))
  allocate(srmatbox3_z_3d(mg_num(1),mg_num(2),Nd))
  allocate(srmatbox4_x_3d(Nd,mg_num(2),mg_num(3)))
  allocate(srmatbox4_y_3d(mg_num(1),Nd,mg_num(3)))
  allocate(srmatbox4_z_3d(mg_num(1),mg_num(2),Nd))

  srmatbox1_x_3d=0.d0
  srmatbox1_y_3d=0.d0
  srmatbox1_z_3d=0.d0
  srmatbox2_x_3d=0.d0
  srmatbox2_y_3d=0.d0
  srmatbox2_z_3d=0.d0
  srmatbox3_x_3d=0.d0
  srmatbox3_y_3d=0.d0
  srmatbox3_z_3d=0.d0
  srmatbox4_x_3d=0.d0
  srmatbox4_y_3d=0.d0
  srmatbox4_z_3d=0.d0

  allocate(scmatbox1_x_3d(Nd,mg_num(2),mg_num(3)))
  allocate(scmatbox1_y_3d(mg_num(1),Nd,mg_num(3)))
  allocate(scmatbox1_z_3d(mg_num(1),mg_num(2),Nd))
  allocate(scmatbox2_x_3d(Nd,mg_num(2),mg_num(3)))
  allocate(scmatbox2_y_3d(mg_num(1),Nd,mg_num(3)))
  allocate(scmatbox2_z_3d(mg_num(1),mg_num(2),Nd))
  allocate(scmatbox3_x_3d(Nd,mg_num(2),mg_num(3)))
  allocate(scmatbox3_y_3d(mg_num(1),Nd,mg_num(3)))
  allocate(scmatbox3_z_3d(mg_num(1),mg_num(2),Nd))
  allocate(scmatbox4_x_3d(Nd,mg_num(2),mg_num(3)))
  allocate(scmatbox4_y_3d(mg_num(1),Nd,mg_num(3)))
  allocate(scmatbox4_z_3d(mg_num(1),mg_num(2),Nd))

  scmatbox1_x_3d=0.d0
  scmatbox1_y_3d=0.d0
  scmatbox1_z_3d=0.d0
  scmatbox2_x_3d=0.d0
  scmatbox2_y_3d=0.d0
  scmatbox2_z_3d=0.d0
  scmatbox3_x_3d=0.d0
  scmatbox3_y_3d=0.d0
  scmatbox3_z_3d=0.d0
  scmatbox4_x_3d=0.d0
  scmatbox4_y_3d=0.d0
  scmatbox4_z_3d=0.d0
  
  if(iSCFRT==1.and.icalcforce==1)then
  
    allocate(srmatbox1_x_5d(Nd,mg_num(2),mg_num(3),1:iobnum,k_sta:k_end))
    allocate(srmatbox1_y_5d(mg_num(1),Nd,mg_num(3),1:iobnum,k_sta:k_end))
    allocate(srmatbox1_z_5d(mg_num(1),mg_num(2),Nd,1:iobnum,k_sta:k_end))
    allocate(srmatbox2_x_5d(Nd,mg_num(2),mg_num(3),1:iobnum,k_sta:k_end))
    allocate(srmatbox2_y_5d(mg_num(1),Nd,mg_num(3),1:iobnum,k_sta:k_end))
    allocate(srmatbox2_z_5d(mg_num(1),mg_num(2),Nd,1:iobnum,k_sta:k_end))
    allocate(srmatbox3_x_5d(Nd,mg_num(2),mg_num(3),1:iobnum,k_sta:k_end))
    allocate(srmatbox3_y_5d(mg_num(1),Nd,mg_num(3),1:iobnum,k_sta:k_end))
    allocate(srmatbox3_z_5d(mg_num(1),mg_num(2),Nd,1:iobnum,k_sta:k_end))
    allocate(srmatbox4_x_5d(Nd,mg_num(2),mg_num(3),1:iobnum,k_sta:k_end))
    allocate(srmatbox4_y_5d(mg_num(1),Nd,mg_num(3),1:iobnum,k_sta:k_end))
    allocate(srmatbox4_z_5d(mg_num(1),mg_num(2),Nd,1:iobnum,k_sta:k_end))
  
    srmatbox1_x_5d=0.d0
    srmatbox1_y_5d=0.d0
    srmatbox1_z_5d=0.d0
    srmatbox2_x_5d=0.d0
    srmatbox2_y_5d=0.d0
    srmatbox2_z_5d=0.d0
    srmatbox3_x_5d=0.d0
    srmatbox3_y_5d=0.d0
    srmatbox3_z_5d=0.d0
    srmatbox4_x_5d=0.d0
    srmatbox4_y_5d=0.d0
    srmatbox4_z_5d=0.d0
    
  else if(iSCFRT==2)then
  
    allocate(scmatbox1_x_5d(Nd,mg_num(2),mg_num(3),1:iobnum,k_sta:k_end))
    allocate(scmatbox1_y_5d(mg_num(1),Nd,mg_num(3),1:iobnum,k_sta:k_end))
    allocate(scmatbox1_z_5d(mg_num(1),mg_num(2),Nd,1:iobnum,k_sta:k_end))
    allocate(scmatbox2_x_5d(Nd,mg_num(2),mg_num(3),1:iobnum,k_sta:k_end))
    allocate(scmatbox2_y_5d(mg_num(1),Nd,mg_num(3),1:iobnum,k_sta:k_end))
    allocate(scmatbox2_z_5d(mg_num(1),mg_num(2),Nd,1:iobnum,k_sta:k_end))
    allocate(scmatbox3_x_5d(Nd,mg_num(2),mg_num(3),1:iobnum,k_sta:k_end))
    allocate(scmatbox3_y_5d(mg_num(1),Nd,mg_num(3),1:iobnum,k_sta:k_end))
    allocate(scmatbox3_z_5d(mg_num(1),mg_num(2),Nd,1:iobnum,k_sta:k_end))
    allocate(scmatbox4_x_5d(Nd,mg_num(2),mg_num(3),1:iobnum,k_sta:k_end))
    allocate(scmatbox4_y_5d(mg_num(1),Nd,mg_num(3),1:iobnum,k_sta:k_end))
    allocate(scmatbox4_z_5d(mg_num(1),mg_num(2),Nd,1:iobnum,k_sta:k_end))
  
    scmatbox1_x_5d=0.d0
    scmatbox1_y_5d=0.d0
    scmatbox1_z_5d=0.d0
    scmatbox2_x_5d=0.d0
    scmatbox2_y_5d=0.d0
    scmatbox2_z_5d=0.d0
    scmatbox3_x_5d=0.d0
    scmatbox3_y_5d=0.d0
    scmatbox3_z_5d=0.d0
    scmatbox4_x_5d=0.d0
    scmatbox4_y_5d=0.d0
    scmatbox4_z_5d=0.d0
  
  end if
  
end subroutine allocate_sendrecv
