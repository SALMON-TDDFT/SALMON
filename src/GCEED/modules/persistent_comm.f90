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
module persistent_comm
  use pack_unpack, only: array_shape

integer,allocatable :: nreqs_rgroupob(:),nreqs_cgroupob(:)

  type(array_shape),public,allocatable :: nshape_groupob(:),nrange_groupob(:,:)

  public :: init_persistent_requests

private
contains
  subroutine init_persistent_requests
    implicit none

    call init_comm_groupob
  end subroutine

  subroutine init_comm_groupob
    use init_sendrecv_sub,    only: iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array
    use salmon_parallel,      only: icomm => nproc_group_orbital
    use salmon_communication, only: comm_send_init, comm_recv_init
    use pack_unpack,          only: create_array_shape
    use scf_data
    implicit none
    integer :: iup,idw,jup,jdw,kup,kdw

    iup=iup_array(1)
    idw=idw_array(1)
    jup=jup_array(1)
    jdw=jdw_array(1)
    kup=kup_array(1)
    kdw=kdw_array(1)

    if(iSCFRT==1.and.icalcforce==1)then
      allocate(nreqs_rgroupob(12))
      nreqs_rgroupob( 1) = comm_send_init(srmatbox1_x_5d,iup,3,icomm)
      nreqs_rgroupob( 2) = comm_recv_init(srmatbox2_x_5d,idw,3,icomm)
      nreqs_rgroupob( 3) = comm_send_init(srmatbox3_x_5d,idw,4,icomm)
      nreqs_rgroupob( 4) = comm_recv_init(srmatbox4_x_5d,iup,4,icomm)
      nreqs_rgroupob( 5) = comm_send_init(srmatbox1_y_5d,jup,5,icomm)
      nreqs_rgroupob( 6) = comm_recv_init(srmatbox2_y_5d,jdw,5,icomm)
      nreqs_rgroupob( 7) = comm_send_init(srmatbox3_y_5d,jdw,6,icomm)
      nreqs_rgroupob( 8) = comm_recv_init(srmatbox4_y_5d,jup,6,icomm)
      nreqs_rgroupob( 9) = comm_send_init(srmatbox1_z_5d,kup,7,icomm)
      nreqs_rgroupob(10) = comm_recv_init(srmatbox2_z_5d,kdw,7,icomm)
      nreqs_rgroupob(11) = comm_send_init(srmatbox3_z_5d,kdw,8,icomm)
      nreqs_rgroupob(12) = comm_recv_init(srmatbox4_z_5d,kup,8,icomm)
    else if(iSCFRT==2.and.nproc_Mxin_mul/=1)then
      allocate(nreqs_cgroupob(12))
      nreqs_cgroupob( 1) = comm_send_init(scmatbox1_x_5d,iup,3,icomm)
      nreqs_cgroupob( 2) = comm_recv_init(scmatbox2_x_5d,idw,3,icomm)
      nreqs_cgroupob( 3) = comm_send_init(scmatbox3_x_5d,idw,4,icomm)
      nreqs_cgroupob( 4) = comm_recv_init(scmatbox4_x_5d,iup,4,icomm)
      nreqs_cgroupob( 5) = comm_send_init(scmatbox1_y_5d,jup,5,icomm)
      nreqs_cgroupob( 6) = comm_recv_init(scmatbox2_y_5d,jdw,5,icomm)
      nreqs_cgroupob( 7) = comm_send_init(scmatbox3_y_5d,jdw,6,icomm)
      nreqs_cgroupob( 8) = comm_recv_init(scmatbox4_y_5d,jup,6,icomm)
      nreqs_cgroupob( 9) = comm_send_init(scmatbox1_z_5d,kup,7,icomm)
      nreqs_cgroupob(10) = comm_recv_init(scmatbox2_z_5d,kdw,7,icomm)
      nreqs_cgroupob(11) = comm_send_init(scmatbox3_z_5d,kdw,8,icomm)
      nreqs_cgroupob(12) = comm_recv_init(scmatbox4_z_5d,kup,8,icomm)
    end if
    allocate(nshape_groupob(5))
    allocate(nrange_groupob(5,12))

    nshape_groupob(1) = create_array_shape(mg_sta(1)-Nd, mg_end(1)+Nd+1)
    nshape_groupob(2) = create_array_shape(mg_sta(2)-Nd, mg_end(2)+Nd)
    nshape_groupob(3) = create_array_shape(mg_sta(3)-Nd, mg_end(3)+Nd)
    nshape_groupob(4) = create_array_shape(1,iobnum)
    nshape_groupob(5) = create_array_shape(1,1)

    nrange_groupob(4,:) = create_array_shape(1,iobnum)
    nrange_groupob(5,:) = create_array_shape(1,1)

    nrange_groupob(1,1) = create_array_shape(mg_end(1)-Nd+1,mg_end(1))
    nrange_groupob(2,1) = create_array_shape(mg_sta(2),     mg_sta(2)+mg_num(2)-1)
    nrange_groupob(3,1) = create_array_shape(mg_sta(3),     mg_sta(3)+mg_num(3)-1)

    nrange_groupob(1,2) = create_array_shape(mg_sta(1),mg_sta(1)+Nd-1)
    nrange_groupob(2,2) = create_array_shape(mg_sta(2),mg_sta(2)+mg_num(2)-1)
    nrange_groupob(3,2) = create_array_shape(mg_sta(3),mg_sta(3)+mg_num(3)-1)

    nrange_groupob(1,3) = create_array_shape(mg_sta(1),     mg_sta(1)+mg_num(1)-1)
    nrange_groupob(2,3) = create_array_shape(mg_end(2)-Nd+1,mg_end(2))
    nrange_groupob(3,3) = create_array_shape(mg_sta(3),     mg_sta(3)+mg_num(3)-1)

    nrange_groupob(1,4) = create_array_shape(mg_sta(1),mg_sta(1)+mg_num(1)-1)
    nrange_groupob(2,4) = create_array_shape(mg_sta(2),mg_sta(2)+Nd-1)
    nrange_groupob(3,4) = create_array_shape(mg_sta(3),mg_sta(3)+mg_num(3)-1)

    nrange_groupob(1,5) = create_array_shape(mg_sta(1),     mg_sta(1)+mg_num(1)-1)
    nrange_groupob(2,5) = create_array_shape(mg_sta(2),     mg_sta(2)+mg_num(2)-1)
    nrange_groupob(3,5) = create_array_shape(mg_end(3)-Nd+1,mg_end(3))

    nrange_groupob(1,6) = create_array_shape(mg_sta(1),mg_sta(1)+mg_num(1)-1)
    nrange_groupob(2,6) = create_array_shape(mg_sta(2),mg_sta(2)+mg_num(2)-1)
    nrange_groupob(3,6) = create_array_shape(mg_sta(3),mg_sta(3)+Nd-1)

    nrange_groupob(1,7) = create_array_shape(mg_sta(1)-Nd,mg_sta(1)-1)
    nrange_groupob(2,7) = create_array_shape(mg_sta(2),   mg_sta(2)+mg_num(2)-1)
    nrange_groupob(3,7) = create_array_shape(mg_sta(3),   mg_sta(3)+mg_num(3)-1)

    nrange_groupob(1,8) = create_array_shape(mg_end(1)+1,mg_end(1)+Nd)
    nrange_groupob(2,8) = create_array_shape(mg_sta(2),  mg_sta(2)+mg_num(2)-1)
    nrange_groupob(3,8) = create_array_shape(mg_sta(3),  mg_sta(3)+mg_num(3)-1)

    nrange_groupob(1,9) = create_array_shape(mg_sta(1),   mg_sta(1)+mg_num(1)-1)
    nrange_groupob(2,9) = create_array_shape(mg_sta(2)-Nd,mg_sta(2)-1)
    nrange_groupob(3,9) = create_array_shape(mg_sta(3),   mg_sta(3)+mg_num(3)-1)

    nrange_groupob(1,10) = create_array_shape(mg_sta(1),  mg_sta(1)+mg_num(1)-1)
    nrange_groupob(2,10) = create_array_shape(mg_end(2)+1,mg_end(2)+Nd)
    nrange_groupob(3,10) = create_array_shape(mg_sta(3),  mg_sta(3)+mg_num(3)-1)

    nrange_groupob(1,11) = create_array_shape(mg_sta(1),   mg_sta(1)+mg_num(1)-1)
    nrange_groupob(2,11) = create_array_shape(mg_sta(2),   mg_sta(2)+mg_num(2)-1)
    nrange_groupob(3,11) = create_array_shape(mg_sta(3)-Nd,mg_sta(3)-1)

    nrange_groupob(1,12) = create_array_shape(mg_sta(1),  mg_sta(1)+mg_num(1)-1)
    nrange_groupob(2,12) = create_array_shape(mg_sta(2),  mg_sta(2)+mg_num(2)-1)
    nrange_groupob(3,12) = create_array_shape(mg_end(3)+1,mg_end(3)+Nd)
  end subroutine
end module
