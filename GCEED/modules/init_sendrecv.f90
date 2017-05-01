! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

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

complex(8),allocatable :: cmatbox1_x_h(:,:,:),cmatbox1_y_h(:,:,:),cmatbox1_z_h(:,:,:)
complex(8),allocatable :: cmatbox2_x_h(:,:,:),cmatbox2_y_h(:,:,:),cmatbox2_z_h(:,:,:)

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE init_updown
!$ use omp_lib
use new_world_sub

implicit none

iup_array(1)=newrank_comm_orbital+1
idw_array(1)=newrank_comm_orbital-1
if(imr(1)==nproc_Mxin(1)-1) iup_array(1)=MPI_PROC_NULL
if(imr(1)==0) idw_array(1)=MPI_PROC_NULL

jup_array(1)=newrank_comm_orbital+nproc_Mxin(1)
jdw_array(1)=newrank_comm_orbital-nproc_Mxin(1)
if(imr(2)==nproc_Mxin(2)-1) jup_array(1)=MPI_PROC_NULL
if(imr(2)==0) jdw_array(1)=MPI_PROC_NULL

kup_array(1)=newrank_comm_orbital+nproc_Mxin(1)*nproc_Mxin(2)
kdw_array(1)=newrank_comm_orbital-nproc_Mxin(1)*nproc_Mxin(2)
if(imr(3)==nproc_Mxin(3)-1) kup_array(1)=MPI_PROC_NULL
if(imr(3)==0) kdw_array(1)=MPI_PROC_NULL

if(isequential==1)then
  if(imr(1)==nproc_Mxin(1)-1.and.imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(2)=MPI_PROC_NULL
  else if(imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(2)=newrank_comm_h+nproc_Mxin_mul_s_dm-nproc_Mxin_s_dm(1)+1
  else
    iup_array(2)=newrank_comm_h+1
  end if
else if(isequential==2)then
  if(imr(1)==nproc_Mxin(1)-1.and.imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(2)=MPI_PROC_NULL
  else if(imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(2)=newrank_comm_h+nproc_Mxin_mul+1-nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  else
    iup_array(2)=newrank_comm_h+nproc_Mxin_mul
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(1)==1.and.nproc_Mxin_s_dm(1)==1)then
    iup_array(4)=MPI_PROC_NULL
  else if(imr(1)==nproc_Mxin(1)-1.and.imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(4)=newrank_comm_h+nproc_Mxin_mul+1-nproc_Mxin_mul*nproc_Mxin_s_dm(1)-nproc_Mxin(1)
  else if(imrs(1)==nproc_Mxin_s_dm(1)-1) then
    iup_array(4)=newrank_comm_h+nproc_Mxin_mul+1-nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  else
    iup_array(4)=newrank_comm_h+nproc_Mxin_mul
  end if
end if

if(isequential==1)then
  if(imr(1)==0.and.imrs(1)==0) then
    idw_array(2)=MPI_PROC_NULL
  else if(imrs(1)==0) then
    idw_array(2)=newrank_comm_h-nproc_Mxin_mul_s_dm+nproc_Mxin_s_dm(1)-1
  else
    idw_array(2)=newrank_comm_h-1
  end if
else if(isequential==2)then
  if(imr(1)==0.and.imrs(1)==0) then
    idw_array(2)=MPI_PROC_NULL
  else if(imrs(1)==0) then
    idw_array(2)=newrank_comm_h-nproc_Mxin_mul-1+nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  else
    idw_array(2)=newrank_comm_h-nproc_Mxin_mul
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(1)==1.and.nproc_Mxin_s_dm(1)==1)then
    idw_array(4)=MPI_PROC_NULL
  else if(imr(1)==0.and.imrs(1)==0) then
    idw_array(4)=newrank_comm_h-nproc_Mxin_mul-1+nproc_Mxin_mul*nproc_Mxin_s_dm(1)+nproc_Mxin(1)
  else if(imrs(1)==0) then
    idw_array(4)=newrank_comm_h-nproc_Mxin_mul-1+nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  else
    idw_array(4)=newrank_comm_h-nproc_Mxin_mul
  end if
end if

if(isequential==1)then
  if(imr(2)==nproc_Mxin(2)-1.and.imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(2)=MPI_PROC_NULL
  else if(imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(2)=newrank_comm_h+nproc_Mxin_mul_s_dm*nproc_Mxin(1)    &
                                    -(nproc_Mxin_s_dm(2)-1)*nproc_Mxin_s_dm(1)
  else
    jup_array(2)=newrank_comm_h+nproc_Mxin_s_dm(1)
  end if
else if(isequential==2)then
  if(imr(2)==nproc_Mxin(2)-1.and.imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(2)=MPI_PROC_NULL
  else if(imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(2)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
            +nproc_Mxin(1)-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    jup_array(2)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(2)==1.and.nproc_Mxin_s_dm(2)==1)then
    jup_array(4)=MPI_PROC_NULL
  else if(imr(2)==nproc_Mxin(2)-1.and.imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(4)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
            +nproc_Mxin(1)-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
            -nproc_Mxin(1)*nproc_Mxin(2)
  else if(imrs(2)==nproc_Mxin_s_dm(2)-1) then
    jup_array(4)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
            +nproc_Mxin(1)-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    jup_array(4)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  end if
end if

if(isequential==1)then
  if(imr(2)==0.and.imrs(2)==0) then
    jdw_array(2)=MPI_PROC_NULL
  else if(imrs(2)==0) then
    jdw_array(2)=newrank_comm_h-nproc_Mxin_mul_s_dm*nproc_Mxin(1)    &
                                    +(nproc_Mxin_s_dm(2)-1)*nproc_Mxin_s_dm(1)
  else
    jdw_array(2)=newrank_comm_h-nproc_Mxin_s_dm(1)
  end if
else if(isequential==2)then
  if(imr(2)==0.and.imrs(2)==0) then
    jdw_array(2)=MPI_PROC_NULL
  else if(imrs(2)==0) then
    jdw_array(2)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
              -nproc_Mxin(1)+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    jdw_array(2)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(2)==1.and.nproc_Mxin_s_dm(2)==1)then
    jdw_array(4)=MPI_PROC_NULL
  else if(imr(2)==0.and.imrs(2)==0) then
    jdw_array(4)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
              -nproc_Mxin(1)+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              +nproc_Mxin(1)*nproc_Mxin(2)
  else if(imrs(2)==0) then
    jdw_array(4)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
              -nproc_Mxin(1)+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    jdw_array(4)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)
  end if
end if


if(isequential==1)then
  if(imr(3)==nproc_Mxin(3)-1.and.imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(2)=MPI_PROC_NULL
  else if(imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(2)=newrank_comm_h+nproc_Mxin_mul_s_dm*nproc_Mxin(1)*nproc_Mxin(2)   &
                                    -(nproc_Mxin_s_dm(3)-1)*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    kup_array(2)=newrank_comm_h+nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
else if(isequential==2)then
  if(imr(3)==nproc_Mxin(3)-1.and.imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(2)=MPI_PROC_NULL
  else if(imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(2)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              +nproc_Mxin(1)*nproc_Mxin(2)   &
              -nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  else
    kup_array(2)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(3)==1.and.nproc_Mxin_s_dm(3)==1)then
    kup_array(4)=MPI_PROC_NULL
  else if(imr(3)==nproc_Mxin(3)-1.and.imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(4)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              +nproc_Mxin(1)*nproc_Mxin(2)   &
              -nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)  &
              -nproc_Mxin_mul
  else if(imrs(3)==nproc_Mxin_s_dm(3)-1) then
    kup_array(4)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              +nproc_Mxin(1)*nproc_Mxin(2)   &
              -nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  else
    kup_array(4)=newrank_comm_h+nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
end if

if(isequential==1)then
  if(imr(3)==0.and.imrs(3)==0) then
    kdw_array(2)=MPI_PROC_NULL
  else if(imrs(3)==0) then
    kdw_array(2)=newrank_comm_h-nproc_Mxin_mul_s_dm*nproc_Mxin(1)*nproc_Mxin(2)   &
                                    +(nproc_Mxin_s_dm(3)-1)*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  else
    kdw_array(2)=newrank_comm_h-nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
else if(isequential==2)then
  if(imr(3)==0.and.imrs(3)==0) then
    kdw_array(2)=MPI_PROC_NULL
  else if(imrs(3)==0) then
    kdw_array(2)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              -nproc_Mxin(1)*nproc_Mxin(2)   &
              +nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  else
    kdw_array(2)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
end if

if(isequential==2)then
  if(nproc_Mxin(3)==1.and.nproc_Mxin_s_dm(3)==1)then
    kdw_array(4)=MPI_PROC_NULL
  else if(imr(3)==0.and.imrs(3)==0) then
    kdw_array(4)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              -nproc_Mxin(1)*nproc_Mxin(2)   &
              +nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)  &
              +nproc_Mxin_mul
  else if(imrs(3)==0) then
    kdw_array(4)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
              -nproc_Mxin(1)*nproc_Mxin(2)   &
              +nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  else
    kdw_array(4)=newrank_comm_h-nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  end if
end if

END SUBROUTINE init_updown

!=====================================================================
subroutine init_itype
!$ use omp_lib

implicit none

integer :: ilap_s
integer :: ibox
integer :: isize1(3),isize2(3)
integer :: isubsize1(3),isubsize2(3)
integer :: istart1(3),istart2(3)
integer :: inum_Mxin2(3,0:nproc-1)

do ilap_s=1,2

  if(ilap_s==1)then
    inum_Mxin2(1:3,0:nproc-1)=inum_Mxin(1:3,0:nproc-1)
  else if(ilap_s==2)then
    inum_Mxin2(1:3,0:nproc-1)=inum_Mxin_s(1:3,0:nproc-1)
  end if

!send from idw to iup

  isize1(1:3)=inum_Mxin2(1:3,myrank)+2*Nd
  isize2(1:3)=inum_Mxin2(1:3,myrank)+2*Nd

  isubsize1(1)=Nd
  isubsize1(2)=inum_Mxin2(2,myrank)
  isubsize1(3)=inum_Mxin2(3,myrank)
  istart1(1)=inum_Mxin2(1,myrank)
  istart1(2)=Nd
  istart1(3)=Nd
  isubsize2(1)=Nd
  isubsize2(2)=inum_Mxin2(2,myrank)
  isubsize2(3)=inum_Mxin2(3,myrank)
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
  istart2(1)=Nd+inum_Mxin2(1,myrank)

  if(ilap_s==1)then
    ibox=3
  else if(ilap_s==2)then
    ibox=15
  end if

!send from jdw to jup

  isize1(1:3)=inum_Mxin2(1:3,myrank)+2*Nd
  isubsize1(1)=inum_Mxin2(1,myrank)
  isubsize1(2)=Nd
  isubsize1(3)=inum_Mxin2(3,myrank)
  istart1(1)=Nd
  istart1(2)=inum_Mxin2(2,myrank)
  istart1(3)=Nd
  isize2(1:3)=inum_Mxin2(1:3,myrank)+2*Nd
  isubsize2(1)=inum_Mxin2(1,myrank)
  isubsize2(2)=Nd
  isubsize2(3)=inum_Mxin2(3,myrank)
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
  istart2(2)=Nd+inum_Mxin2(2,myrank)

  if(ilap_s==1)then
    ibox=7
  else if(ilap_s==2)then
    ibox=19
  end if

!send from kup to kdw

  isubsize1(1)=inum_Mxin2(1,myrank)
  isubsize1(2)=inum_Mxin2(2,myrank)
  isubsize1(3)=Nd
  istart1(1)=Nd
  istart1(2)=Nd
  istart1(3)=inum_Mxin2(3,myrank)
  isubsize2(1)=inum_Mxin2(1,myrank)
  isubsize2(2)=inum_Mxin2(2,myrank)
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
  istart2(3)=Nd+inum_Mxin2(3,myrank)

  if(ilap_s==1)then
    ibox=11
  else if(ilap_s==2)then
    ibox=23
  end if

end do

end subroutine init_itype

!======================================================================
subroutine init_sendrecv_matrix
integer :: inum_Mxin2(3,0:nproc-1)

inum_Mxin2(1:3,0:nproc-1)=inum_Mxin_s(1:3,0:nproc-1)

allocate(rmatbox1_x_s(Nd,inum_Mxin2(2,myrank),inum_Mxin2(3,myrank)))
allocate(rmatbox1_y_s(inum_Mxin2(1,myrank),Nd,inum_Mxin2(3,myrank)))
allocate(rmatbox1_z_s(inum_Mxin2(1,myrank),inum_Mxin2(2,myrank),Nd))
allocate(rmatbox2_x_s(Nd,inum_Mxin2(2,myrank),inum_Mxin2(3,myrank)))
allocate(rmatbox2_y_s(inum_Mxin2(1,myrank),Nd,inum_Mxin2(3,myrank)))
allocate(rmatbox2_z_s(inum_Mxin2(1,myrank),inum_Mxin2(2,myrank),Nd))

rmatbox1_x_s=0.d0
rmatbox1_y_s=0.d0
rmatbox1_z_s=0.d0
rmatbox2_x_s=0.d0
rmatbox2_y_s=0.d0
rmatbox2_z_s=0.d0

inum_Mxin2(1:3,0:nproc-1)=inum_Mxin(1:3,0:nproc-1)

allocate(cmatbox1_x_s(Nd,inum_Mxin2(2,myrank),inum_Mxin2(3,myrank)))
allocate(cmatbox1_y_s(inum_Mxin2(1,myrank),Nd,inum_Mxin2(3,myrank)))
allocate(cmatbox1_z_s(inum_Mxin2(1,myrank),inum_Mxin2(2,myrank),Nd))
allocate(cmatbox2_x_s(Nd,inum_Mxin2(2,myrank),inum_Mxin2(3,myrank)))
allocate(cmatbox2_y_s(inum_Mxin2(1,myrank),Nd,inum_Mxin2(3,myrank)))
allocate(cmatbox2_z_s(inum_Mxin2(1,myrank),inum_Mxin2(2,myrank),Nd))

cmatbox1_x_s=0.d0
cmatbox1_y_s=0.d0
cmatbox1_z_s=0.d0
cmatbox2_x_s=0.d0
cmatbox2_y_s=0.d0
cmatbox2_z_s=0.d0

inum_Mxin2(1:3,0:nproc-1)=inum_Mxin_s(1:3,0:nproc-1)

allocate(rmatbox1_x_h(Ndh,inum_Mxin2(2,myrank),inum_Mxin2(3,myrank)))
allocate(rmatbox1_y_h(inum_Mxin2(1,myrank),Ndh,inum_Mxin2(3,myrank)))
allocate(rmatbox1_z_h(inum_Mxin2(1,myrank),inum_Mxin2(2,myrank),Ndh))
allocate(rmatbox2_x_h(Ndh,inum_Mxin2(2,myrank),inum_Mxin2(3,myrank)))
allocate(rmatbox2_y_h(inum_Mxin2(1,myrank),Ndh,inum_Mxin2(3,myrank)))
allocate(rmatbox2_z_h(inum_Mxin2(1,myrank),inum_Mxin2(2,myrank),Ndh))

rmatbox1_x_h=0.d0
rmatbox1_y_h=0.d0
rmatbox1_z_h=0.d0
rmatbox2_x_h=0.d0
rmatbox2_y_h=0.d0
rmatbox2_z_h=0.d0

allocate(cmatbox1_x_h(Ndh,inum_Mxin2(2,myrank),inum_Mxin2(3,myrank)))
allocate(cmatbox1_y_h(inum_Mxin2(1,myrank),Ndh,inum_Mxin2(3,myrank)))
allocate(cmatbox1_z_h(inum_Mxin2(1,myrank),inum_Mxin2(2,myrank),Ndh))
allocate(cmatbox2_x_h(Ndh,inum_Mxin2(2,myrank),inum_Mxin2(3,myrank)))
allocate(cmatbox2_y_h(inum_Mxin2(1,myrank),Ndh,inum_Mxin2(3,myrank)))
allocate(cmatbox2_z_h(inum_Mxin2(1,myrank),inum_Mxin2(2,myrank),Ndh))

cmatbox1_x_h=0.d0
cmatbox1_y_h=0.d0
cmatbox1_z_h=0.d0
cmatbox2_x_h=0.d0
cmatbox2_y_h=0.d0
cmatbox2_z_h=0.d0


end subroutine init_sendrecv_matrix

END MODULE init_sendrecv_sub

