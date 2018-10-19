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

  integer,public,allocatable :: nreqs_rh(:,:)
  integer,public,allocatable :: nreqs_rorbital(:),nreqs_corbital(:)
  integer,public,allocatable :: nreqs_rgroupob(:),nreqs_cgroupob(:)

  type(array_shape),public,allocatable :: nshape_orbital(:)
  type(array_shape),public,allocatable :: nshape_groupob(:),nrange_groupob(:,:)

  public :: init_persistent_requests

private
contains
  subroutine init_persistent_requests
    implicit none

    call init_comm_korbital
    call init_comm_groupob
    call init_comm_h
  end subroutine

  subroutine init_comm_korbital
    use init_sendrecv_sub,    only: iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array
    use salmon_parallel,      only: icomm => nproc_group_korbital
    use salmon_communication, only: comm_send_init, comm_recv_init, comm_proc_null
    use pack_unpack,          only: create_array_shape
    use scf_data
    implicit none
    integer :: iup,idw,jup,jdw,kup,kdw

#ifdef SALMON_USE_MPI
    iup=iup_array(1)
    idw=idw_array(1)
    jup=jup_array(1)
    jdw=jdw_array(1)
    kup=kup_array(1)
    kdw=kdw_array(1)
#else
    iup=comm_proc_null
    idw=comm_proc_null
    jup=comm_proc_null
    jdw=comm_proc_null
    kup=comm_proc_null
    kdw=comm_proc_null
#endif

    allocate(nreqs_rorbital(12))
    nreqs_rorbital( 1) = comm_send_init(srmatbox1_x_3d,iup,3,icomm)
    nreqs_rorbital( 2) = comm_recv_init(srmatbox2_x_3d,idw,3,icomm)
    nreqs_rorbital( 3) = comm_send_init(srmatbox3_x_3d,idw,4,icomm)
    nreqs_rorbital( 4) = comm_recv_init(srmatbox4_x_3d,iup,4,icomm)
    nreqs_rorbital( 5) = comm_send_init(srmatbox1_y_3d,jup,5,icomm)
    nreqs_rorbital( 6) = comm_recv_init(srmatbox2_y_3d,jdw,5,icomm)
    nreqs_rorbital( 7) = comm_send_init(srmatbox3_y_3d,jdw,6,icomm)
    nreqs_rorbital( 8) = comm_recv_init(srmatbox4_y_3d,jup,6,icomm)
    nreqs_rorbital( 9) = comm_send_init(srmatbox1_z_3d,kup,7,icomm)
    nreqs_rorbital(10) = comm_recv_init(srmatbox2_z_3d,kdw,7,icomm)
    nreqs_rorbital(11) = comm_send_init(srmatbox3_z_3d,kdw,8,icomm)
    nreqs_rorbital(12) = comm_recv_init(srmatbox4_z_3d,kup,8,icomm)

    allocate(nreqs_corbital(12))
    nreqs_corbital( 1) = comm_send_init(scmatbox1_x_3d,iup,3,icomm)
    nreqs_corbital( 2) = comm_recv_init(scmatbox2_x_3d,idw,3,icomm)
    nreqs_corbital( 3) = comm_send_init(scmatbox3_x_3d,idw,4,icomm)
    nreqs_corbital( 4) = comm_recv_init(scmatbox4_x_3d,iup,4,icomm)
    nreqs_corbital( 5) = comm_send_init(scmatbox1_y_3d,jup,5,icomm)
    nreqs_corbital( 6) = comm_recv_init(scmatbox2_y_3d,jdw,5,icomm)
    nreqs_corbital( 7) = comm_send_init(scmatbox3_y_3d,jdw,6,icomm)
    nreqs_corbital( 8) = comm_recv_init(scmatbox4_y_3d,jup,6,icomm)
    nreqs_corbital( 9) = comm_send_init(scmatbox1_z_3d,kup,7,icomm)
    nreqs_corbital(10) = comm_recv_init(scmatbox2_z_3d,kdw,7,icomm)
    nreqs_corbital(11) = comm_send_init(scmatbox3_z_3d,kdw,8,icomm)
    nreqs_corbital(12) = comm_recv_init(scmatbox4_z_3d,kup,8,icomm)

    allocate(nshape_orbital(3))
    nshape_orbital(1) = create_array_shape(mg_sta(1)-Nd, mg_end(1)+Nd)
    nshape_orbital(2) = create_array_shape(mg_sta(2)-Nd, mg_end(2)+Nd)
    nshape_orbital(3) = create_array_shape(mg_sta(3)-Nd, mg_end(3)+Nd)
  end subroutine

  subroutine init_comm_groupob
    use init_sendrecv_sub,    only: iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array
    use salmon_parallel,      only: icomm => nproc_group_korbital
    use salmon_communication, only: comm_send_init, comm_recv_init, comm_proc_null
    use pack_unpack,          only: create_array_shape
    use scf_data
    implicit none
    integer :: iup,idw,jup,jdw,kup,kdw

#ifdef SALMON_USE_MPI
    iup=iup_array(1)
    idw=idw_array(1)
    jup=jup_array(1)
    jdw=jdw_array(1)
    kup=kup_array(1)
    kdw=kdw_array(1)
#else
    iup=comm_proc_null
    idw=comm_proc_null
    jup=comm_proc_null
    jdw=comm_proc_null
    kup=comm_proc_null
    kdw=comm_proc_null
#endif

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
    nshape_groupob(5) = create_array_shape(k_sta,k_end)

    nrange_groupob(4,:) = create_array_shape(1,iobnum)
    nrange_groupob(5,:) = create_array_shape(k_sta,k_end)

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

  subroutine init_comm_h
    implicit none

    allocate(nreqs_rh(12,4))
    call init_reqs_h(nreqs_rh(:,1),1)
    call init_reqs_h(nreqs_rh(:,2),2)
    call init_reqs_h(nreqs_rh(:,4),4)
  end subroutine

  subroutine init_reqs_h(ireqs,itype)
    use init_sendrecv_sub
    use salmon_parallel,      only: nproc_group_global, nproc_group_h
    use salmon_communication, only: comm_send_init, comm_recv_init, comm_proc_null
    implicit none
    integer, intent(out) :: ireqs(12)
    integer, intent(in)  :: itype
    integer :: icomm
    integer :: iup,idw,jup,jdw,kup,kdw

#ifdef SALMON_USE_MPI
    iup=iup_array(itype)
    idw=idw_array(itype)
    jup=jup_array(itype)
    jdw=jdw_array(itype)
    kup=kup_array(itype)
    kdw=kdw_array(itype)
#else
    iup=comm_proc_null
    idw=comm_proc_null
    jup=comm_proc_null
    jdw=comm_proc_null
    kup=comm_proc_null
    kdw=comm_proc_null
#endif

    if (itype == 1) then
      icomm = nproc_group_global
    else
      icomm = nproc_group_h
    end if

    ireqs( 1) = comm_send_init(rmatbox1_x_h,iup,3,icomm)
    ireqs( 2) = comm_recv_init(rmatbox2_x_h,idw,3,icomm)
    ireqs( 3) = comm_send_init(rmatbox3_x_h,idw,4,icomm)
    ireqs( 4) = comm_recv_init(rmatbox4_x_h,iup,4,icomm)
    ireqs( 5) = comm_send_init(rmatbox1_y_h,jup,5,icomm)
    ireqs( 6) = comm_recv_init(rmatbox2_y_h,jdw,5,icomm)
    ireqs( 7) = comm_send_init(rmatbox3_y_h,jdw,6,icomm)
    ireqs( 8) = comm_recv_init(rmatbox4_y_h,jup,6,icomm)
    ireqs( 9) = comm_send_init(rmatbox1_z_h,kup,7,icomm)
    ireqs(10) = comm_recv_init(rmatbox2_z_h,kdw,7,icomm)
    ireqs(11) = comm_send_init(rmatbox3_z_h,kdw,8,icomm)
    ireqs(12) = comm_recv_init(rmatbox4_z_h,kup,8,icomm)
  end subroutine
end module
