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
MODULE init_sendrecv_sub

use scf_data
implicit none 
integer :: iup_array(4)
integer :: idw_array(4)
integer :: jup_array(4)
integer :: jdw_array(4)
integer :: kup_array(4)
integer :: kdw_array(4)

integer :: itype_r(24)
integer :: itype_c(24)

real(8),allocatable :: rmatbox1_x_s(:,:,:),rmatbox1_y_s(:,:,:),rmatbox1_z_s(:,:,:)
real(8),allocatable :: rmatbox2_x_s(:,:,:),rmatbox2_y_s(:,:,:),rmatbox2_z_s(:,:,:)

complex(8),allocatable :: cmatbox1_x_s(:,:,:),cmatbox1_y_s(:,:,:),cmatbox1_z_s(:,:,:)
complex(8),allocatable :: cmatbox2_x_s(:,:,:),cmatbox2_y_s(:,:,:),cmatbox2_z_s(:,:,:)

real(8),allocatable :: rmatbox1_x_h(:,:,:),rmatbox1_y_h(:,:,:),rmatbox1_z_h(:,:,:)
real(8),allocatable :: rmatbox2_x_h(:,:,:),rmatbox2_y_h(:,:,:),rmatbox2_z_h(:,:,:)
real(8),allocatable :: rmatbox3_x_h(:,:,:),rmatbox3_y_h(:,:,:),rmatbox3_z_h(:,:,:)
real(8),allocatable :: rmatbox4_x_h(:,:,:),rmatbox4_y_h(:,:,:),rmatbox4_z_h(:,:,:)

complex(8),allocatable :: cmatbox1_x_h(:,:,:),cmatbox1_y_h(:,:,:),cmatbox1_z_h(:,:,:)
complex(8),allocatable :: cmatbox2_x_h(:,:,:),cmatbox2_y_h(:,:,:),cmatbox2_z_h(:,:,:)

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE init_updown
use inputoutput, only: iperiodic
use salmon_parallel
use salmon_communication, only: comm_proc_null
use new_world_sub

implicit none

iup_array(1)=nproc_id_korbital+1
idw_array(1)=nproc_id_korbital-1
select case(iperiodic)
case(0)
  if(imr(1)==nproc_Mxin(1)-1) iup_array(1)=comm_proc_null
  if(imr(1)==0) idw_array(1)=comm_proc_null
case(3)
  if(imr(1)==nproc_Mxin(1)-1) then
    iup_array(1)=nproc_id_korbital-(nproc_Mxin(1)-1)
  end if
  if(imr(1)==0) then
    idw_array(1)=nproc_id_korbital+(nproc_Mxin(1)-1)
  end if
end select

jup_array(1)=nproc_id_korbital+nproc_Mxin(1)
jdw_array(1)=nproc_id_korbital-nproc_Mxin(1)
select case(iperiodic)
case(0)
  if(imr(2)==nproc_Mxin(2)-1) jup_array(1)=comm_proc_null
  if(imr(2)==0) jdw_array(1)=comm_proc_null
case(3)
  if(imr(2)==nproc_Mxin(2)-1) then
    jup_array(1)=nproc_id_korbital-(nproc_Mxin(2)-1)*nproc_Mxin(1)
  end if
  if(imr(2)==0) then
    jdw_array(1)=nproc_id_korbital+(nproc_Mxin(2)-1)*nproc_Mxin(1)
  end if
end select

kup_array(1)=nproc_id_korbital+nproc_Mxin(1)*nproc_Mxin(2)
kdw_array(1)=nproc_id_korbital-nproc_Mxin(1)*nproc_Mxin(2)
select case(iperiodic)
case(0)
  if(imr(3)==nproc_Mxin(3)-1) kup_array(1)=comm_proc_null
  if(imr(3)==0) kdw_array(1)=comm_proc_null
case(3)
  if(imr(3)==nproc_Mxin(3)-1) then
    kup_array(1)=nproc_id_korbital-(nproc_Mxin(3)-1)*nproc_Mxin(1)*nproc_Mxin(2)
  end if
  if(imr(3)==0) then
    kdw_array(1)=nproc_id_korbital+(nproc_Mxin(3)-1)*nproc_Mxin(1)*nproc_Mxin(2)
  end if
end select

if(isequential==1)then
  if(imr(1)==nproc_Mxin(1)-1.and.imrs(1)==nproc_Mxin_s_dm(1)-1) then
    select case(iperiodic)
    case(0)
      iup_array(2)=comm_proc_null
    case(3)
      iup_array(2)=nproc_id_h+nproc_Mxin_mul_s_dm-nproc_Mxin_s_dm(1)+1-nproc_Mxin_mul_s_dm*nproc_Mxin(1)
    end select
  else if(imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(2)=nproc_id_h+nproc_Mxin_mul_s_dm-nproc_Mxin_s_dm(1)+1
  else
    iup_array(2)=nproc_id_h+1
  end if
else if(isequential==2)then
  if(imr(1)==nproc_Mxin(1)-1.and.imrs(1)==nproc_Mxin_s_dm(1)-1) then
    select case(iperiodic)
    case(0)
      iup_array(2)=comm_proc_null
    case(3)
      iup_array(2)=nproc_id_h+nproc_Mxin_mul+1-nproc_Mxin_mul*nproc_Mxin_s_dm(1)-nproc_Mxin(1)
    end select
  else if(imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(2)=nproc_id_h+nproc_Mxin_mul+1-nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  else
    iup_array(2)=nproc_id_h+nproc_Mxin_mul
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(1)==1.and.nproc_Mxin_s_dm(1)==1)then
    iup_array(4)=comm_proc_null
  else if(imr(1)==nproc_Mxin(1)-1.and.imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(4)=nproc_id_h+nproc_Mxin_mul+1-nproc_Mxin_mul*nproc_Mxin_s_dm(1)-nproc_Mxin(1)
  else if(imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(4)=nproc_id_h+nproc_Mxin_mul+1-nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  else
    iup_array(4)=nproc_id_h+nproc_Mxin_mul
  end if
end if

if(isequential==1)then
  if(imr(1)==0.and.imrs(1)==0) then
    select case(iperiodic)
    case(0)
      idw_array(2)=comm_proc_null
    case(3)
      idw_array(2)=nproc_id_h-nproc_Mxin_mul-1+nproc_Mxin_mul*nproc_Mxin_s_dm(1)+nproc_Mxin(1)
    end select
  else if(imrs(1)==0) then
    idw_array(2)=nproc_id_h-nproc_Mxin_mul_s_dm+nproc_Mxin_s_dm(1)-1
  else
    idw_array(2)=nproc_id_h-1
  end if
else if(isequential==2)then
  if(imr(1)==0.and.imrs(1)==0) then
    idw_array(2)=comm_proc_null
  else if(imrs(1)==0) then
    idw_array(2)=nproc_id_h-nproc_Mxin_mul-1+nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  else
    idw_array(2)=nproc_id_h-nproc_Mxin_mul
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(1)==1.and.nproc_Mxin_s_dm(1)==1)then
    idw_array(4)=comm_proc_null
  else if(imr(1)==0.and.imrs(1)==0) then
    idw_array(4)=nproc_id_h-nproc_Mxin_mul-1+nproc_Mxin_mul*nproc_Mxin_s_dm(1)+nproc_Mxin(1)
  else if(imrs(1)==0) then
    idw_array(4)=nproc_id_h-nproc_Mxin_mul-1+nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  else
    idw_array(4)=nproc_id_h-nproc_Mxin_mul
  end if
end if

if(isequential==1)then
  if(imr(2)==nproc_Mxin(2)-1.and.imrs(2)==nproc_Mxin_s_dm(2)-1) then
    select case(iperiodic)
    case(0)
      jup_array(2)=comm_proc_null
    case(3)
      jup_array(2)=nproc_id_h+nproc_Mxin_mul_s_dm*nproc_Mxin(1)    &
                                    -(nproc_Mxin_s_dm(2)-1)*nproc_Mxin_s_dm(1)-nproc_Mxin_mul_s_dm*nproc_Mxin(1)*nproc_Mxin(2)
    end select 
  else if(imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(2)=nproc_id_h+nproc_Mxin_mul_s_dm*nproc_Mxin(1)    &
                                    -(nproc_Mxin_s_dm(2)-1)*nproc_Mxin_s_dm(1)
  else
    jup_array(2)=nproc_id_h+nproc_Mxin_s_dm(1)
  end if
else if(isequential==2)then
  if(imr(2)==nproc_Mxin(2)-1.and.imrs(2)==nproc_Mxin_s_dm(2)-1) then
    select case(iperiodic)
    case(0)
      jup_array(2)=comm_proc_null
    case(3)
      jup_array(2)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
            +nproc_Mxin(1)-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)-nproc_Mxin(1)*nproc_Mxin(2)
    end select
  else if(imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(2)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
            +nproc_Mxin(1)-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    jup_array(2)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(2)==1.and.nproc_Mxin_s_dm(2)==1)then
    jup_array(4)=comm_proc_null
  else if(imr(2)==nproc_Mxin(2)-1.and.imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(4)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
            +nproc_Mxin(1)-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
            -nproc_Mxin(1)*nproc_Mxin(2)
  else if(imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(4)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
            +nproc_Mxin(1)-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    jup_array(4)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  end if
end if

if(isequential==1)then
  if(imr(2)==0.and.imrs(2)==0) then
    select case(iperiodic)
    case(0)
      jdw_array(2)=comm_proc_null
    case(3)
      jdw_array(2)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
              -nproc_Mxin(1)+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)+nproc_Mxin(1)*nproc_Mxin(2)
    end select
  else if(imrs(2)==0) then
    jdw_array(2)=nproc_id_h-nproc_Mxin_mul_s_dm*nproc_Mxin(1)    &
                                    +(nproc_Mxin_s_dm(2)-1)*nproc_Mxin_s_dm(1)
  else
    jdw_array(2)=nproc_id_h-nproc_Mxin_s_dm(1)
  end if
else if(isequential==2)then
  if(imr(2)==0.and.imrs(2)==0) then
    jdw_array(2)=comm_proc_null
  else if(imrs(2)==0) then
    jdw_array(2)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
              -nproc_Mxin(1)+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    jdw_array(2)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(2)==1.and.nproc_Mxin_s_dm(2)==1)then
    jdw_array(4)=comm_proc_null
  else if(imr(2)==0.and.imrs(2)==0) then
    jdw_array(4)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
              -nproc_Mxin(1)+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              +nproc_Mxin(1)*nproc_Mxin(2)
  else if(imrs(2)==0) then
    jdw_array(4)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
              -nproc_Mxin(1)+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    jdw_array(4)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  end if
end if


if(isequential==1)then
  if(imr(3)==nproc_Mxin(3)-1.and.imrs(3)==nproc_Mxin_s_dm(3)-1) then
    select case(iperiodic)
    case(0)
      kup_array(2)=comm_proc_null
    case(3)
      kup_array(2)=nproc_id_h+nproc_Mxin_mul_s_dm*nproc_Mxin(1)*nproc_Mxin(2) &
                                    -(nproc_Mxin_s_dm(3)-1)*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2) &
                                    -nproc_Mxin_mul_s_dm*nproc_Mxin_mul
    end select
  else if(imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(2)=nproc_id_h+nproc_Mxin_mul_s_dm*nproc_Mxin(1)*nproc_Mxin(2)   &
                                    -(nproc_Mxin_s_dm(3)-1)*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    kup_array(2)=nproc_id_h+nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
else if(isequential==2)then
  if(imr(3)==nproc_Mxin(3)-1.and.imrs(3)==nproc_Mxin_s_dm(3)-1) then
    select case(iperiodic)
    case(0)
      kup_array(2)=comm_proc_null
    case(3)
      kup_array(2)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2) &
              +nproc_Mxin(1)*nproc_Mxin(2)   &
              -nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)-nproc_Mxin_mul
    end select
  else if(imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(2)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              +nproc_Mxin(1)*nproc_Mxin(2)   &
              -nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  else
    kup_array(2)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(3)==1.and.nproc_Mxin_s_dm(3)==1)then
    kup_array(4)=comm_proc_null
  else if(imr(3)==nproc_Mxin(3)-1.and.imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(4)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              +nproc_Mxin(1)*nproc_Mxin(2)   &
              -nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)  &
              -nproc_Mxin_mul
  else if(imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(4)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              +nproc_Mxin(1)*nproc_Mxin(2)   &
              -nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  else
    kup_array(4)=nproc_id_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
end if

if(isequential==1)then
  if(imr(3)==0.and.imrs(3)==0) then
    select case(iperiodic)
    case(0)
      kdw_array(2)=comm_proc_null
    case(3)
      kdw_array(2)=nproc_id_h-nproc_Mxin_mul_s_dm*nproc_Mxin(1)*nproc_Mxin(2) &
                                    +(nproc_Mxin_s_dm(3)-1)*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2) &
                                    +nproc_Mxin_mul_s_dm*nproc_Mxin_mul
    end select
  else if(imrs(3)==0) then
    kdw_array(2)=nproc_id_h-nproc_Mxin_mul_s_dm*nproc_Mxin(1)*nproc_Mxin(2)   &
                                    +(nproc_Mxin_s_dm(3)-1)*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    kdw_array(2)=nproc_id_h-nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
else if(isequential==2)then
  if(imr(3)==0.and.imrs(3)==0) then
    select case(iperiodic)
    case(0)
      kdw_array(2)=comm_proc_null
    case(3)
      kdw_array(2)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2) &
              -nproc_Mxin(1)*nproc_Mxin(2)   &
              +nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)+nproc_Mxin_mul
    end select
  else if(imrs(3)==0) then
    kdw_array(2)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              -nproc_Mxin(1)*nproc_Mxin(2)   &
              +nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  else
    kdw_array(2)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(3)==1.and.nproc_Mxin_s_dm(3)==1)then
    kdw_array(4)=comm_proc_null
  else if(imr(3)==0.and.imrs(3)==0) then
    kdw_array(4)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              -nproc_Mxin(1)*nproc_Mxin(2)   &
              +nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)  &
              +nproc_Mxin_mul
  else if(imrs(3)==0) then
    kdw_array(4)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              -nproc_Mxin(1)*nproc_Mxin(2)   &
              +nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  else
    kdw_array(4)=nproc_id_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
end if

END SUBROUTINE init_updown

!=====================================================================
subroutine init_itype
use salmon_parallel, only: nproc_size_global, nproc_id_global

implicit none

integer :: ilap_s
integer :: ibox
integer :: isize1(3),isize2(3)
integer :: isubsize1(3),isubsize2(3)
integer :: istart1(3),istart2(3)
integer :: inum_Mxin2(3,0:nproc_size_global-1)

do ilap_s=1,2

  if(ilap_s==1)then
    inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin(1:3,0:nproc_size_global-1)
  else if(ilap_s==2)then
    inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin_s(1:3,0:nproc_size_global-1)
  end if

!send from idw to iup

  isize1(1:3)=inum_Mxin2(1:3,nproc_id_global)+2*Nd
  isize2(1:3)=inum_Mxin2(1:3,nproc_id_global)+2*Nd

  isubsize1(1)=Nd
  isubsize1(2)=inum_Mxin2(2,nproc_id_global)
  isubsize1(3)=inum_Mxin2(3,nproc_id_global)
  istart1(1)=inum_Mxin2(1,nproc_id_global)
  istart1(2)=Nd
  istart1(3)=Nd
  isubsize2(1)=Nd
  isubsize2(2)=inum_Mxin2(2,nproc_id_global)
  isubsize2(3)=inum_Mxin2(3,nproc_id_global)
  istart2(1)=0
  istart2(2)=Nd
  istart2(3)=Nd

  if(ilap_s==1)then
    ibox=1
  else if(ilap_s==2)then
    ibox=13
  end if

!send from iup to idw

  istart1(1)=Nd
  istart2(1)=Nd+inum_Mxin2(1,nproc_id_global)

  if(ilap_s==1)then
    ibox=3
  else if(ilap_s==2)then
    ibox=15
  end if

!send from jdw to jup

  isize1(1:3)=inum_Mxin2(1:3,nproc_id_global)+2*Nd
  isubsize1(1)=inum_Mxin2(1,nproc_id_global)
  isubsize1(2)=Nd
  isubsize1(3)=inum_Mxin2(3,nproc_id_global)
  istart1(1)=Nd
  istart1(2)=inum_Mxin2(2,nproc_id_global)
  istart1(3)=Nd
  isize2(1:3)=inum_Mxin2(1:3,nproc_id_global)+2*Nd
  isubsize2(1)=inum_Mxin2(1,nproc_id_global)
  isubsize2(2)=Nd
  isubsize2(3)=inum_Mxin2(3,nproc_id_global)
  istart2(1)=Nd
  istart2(2)=0
  istart2(3)=Nd

  if(ilap_s==1)then
    ibox=5
  else if(ilap_s==2)then
    ibox=17
  end if

!send from jup to jdw
  
  istart1(2)=Nd
  istart2(2)=Nd+inum_Mxin2(2,nproc_id_global)

  if(ilap_s==1)then
    ibox=7
  else if(ilap_s==2)then
    ibox=19
  end if

!send from kup to kdw

  isubsize1(1)=inum_Mxin2(1,nproc_id_global)
  isubsize1(2)=inum_Mxin2(2,nproc_id_global)
  isubsize1(3)=Nd
  istart1(1)=Nd
  istart1(2)=Nd
  istart1(3)=inum_Mxin2(3,nproc_id_global)
  isubsize2(1)=inum_Mxin2(1,nproc_id_global)
  isubsize2(2)=inum_Mxin2(2,nproc_id_global)
  isubsize2(3)=Nd
  istart2(1)=Nd
  istart2(2)=Nd
  istart2(3)=0

  if(ilap_s==1)then
    ibox=9
  else if(ilap_s==2)then
    ibox=21
  end if

!send from kup to kdw

  istart1(3)=Nd
  istart2(3)=Nd+inum_Mxin2(3,nproc_id_global)

  if(ilap_s==1)then
    ibox=11
  else if(ilap_s==2)then
    ibox=23
  end if

end do

end subroutine init_itype

!======================================================================
subroutine init_sendrecv_matrix
use salmon_parallel, only: nproc_size_global, nproc_id_global
integer :: inum_Mxin2(3,0:nproc_size_global-1)

inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin_s(1:3,0:nproc_size_global-1)

allocate(rmatbox1_x_s(Nd,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox1_y_s(inum_Mxin2(1,nproc_id_global),Nd,inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox1_z_s(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Nd))
allocate(rmatbox2_x_s(Nd,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox2_y_s(inum_Mxin2(1,nproc_id_global),Nd,inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox2_z_s(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Nd))

rmatbox1_x_s=0.d0
rmatbox1_y_s=0.d0
rmatbox1_z_s=0.d0
rmatbox2_x_s=0.d0
rmatbox2_y_s=0.d0
rmatbox2_z_s=0.d0

inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin(1:3,0:nproc_size_global-1)

allocate(cmatbox1_x_s(Nd,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(cmatbox1_y_s(inum_Mxin2(1,nproc_id_global),Nd,inum_Mxin2(3,nproc_id_global)))
allocate(cmatbox1_z_s(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Nd))
allocate(cmatbox2_x_s(Nd,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(cmatbox2_y_s(inum_Mxin2(1,nproc_id_global),Nd,inum_Mxin2(3,nproc_id_global)))
allocate(cmatbox2_z_s(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Nd))

cmatbox1_x_s=0.d0
cmatbox1_y_s=0.d0
cmatbox1_z_s=0.d0
cmatbox2_x_s=0.d0
cmatbox2_y_s=0.d0
cmatbox2_z_s=0.d0

inum_Mxin2(1:3,0:nproc_size_global-1)=inum_Mxin_s(1:3,0:nproc_size_global-1)

allocate(rmatbox1_x_h(Ndh,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox1_y_h(inum_Mxin2(1,nproc_id_global),Ndh,inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox1_z_h(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Ndh))
allocate(rmatbox2_x_h(Ndh,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox2_y_h(inum_Mxin2(1,nproc_id_global),Ndh,inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox2_z_h(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Ndh))
allocate(rmatbox3_x_h(Ndh,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox3_y_h(inum_Mxin2(1,nproc_id_global),Ndh,inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox3_z_h(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Ndh))
allocate(rmatbox4_x_h(Ndh,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox4_y_h(inum_Mxin2(1,nproc_id_global),Ndh,inum_Mxin2(3,nproc_id_global)))
allocate(rmatbox4_z_h(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Ndh))

rmatbox1_x_h=0.d0
rmatbox1_y_h=0.d0
rmatbox1_z_h=0.d0
rmatbox2_x_h=0.d0
rmatbox2_y_h=0.d0
rmatbox2_z_h=0.d0
rmatbox3_x_h=0.d0
rmatbox3_y_h=0.d0
rmatbox3_z_h=0.d0
rmatbox4_x_h=0.d0
rmatbox4_y_h=0.d0
rmatbox4_z_h=0.d0

allocate(cmatbox1_x_h(Ndh,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(cmatbox1_y_h(inum_Mxin2(1,nproc_id_global),Ndh,inum_Mxin2(3,nproc_id_global)))
allocate(cmatbox1_z_h(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Ndh))
allocate(cmatbox2_x_h(Ndh,inum_Mxin2(2,nproc_id_global),inum_Mxin2(3,nproc_id_global)))
allocate(cmatbox2_y_h(inum_Mxin2(1,nproc_id_global),Ndh,inum_Mxin2(3,nproc_id_global)))
allocate(cmatbox2_z_h(inum_Mxin2(1,nproc_id_global),inum_Mxin2(2,nproc_id_global),Ndh))

cmatbox1_x_h=0.d0
cmatbox1_y_h=0.d0
cmatbox1_z_h=0.d0
cmatbox2_x_h=0.d0
cmatbox2_y_h=0.d0
cmatbox2_z_h=0.d0


end subroutine init_sendrecv_matrix

END MODULE init_sendrecv_sub

