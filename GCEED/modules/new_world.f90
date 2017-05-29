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
MODULE new_world_sub

use scf_data

implicit none
integer :: newprocs_comm_spin,newrank_comm_spin,newworld_comm_spin
integer :: newprocs_comm_grid,newrank_comm_grid,newworld_comm_grid
integer :: newprocs_comm_orbital,newrank_comm_orbital,newworld_comm_orbital
integer :: newprocs_comm_h,newrank_comm_h,newworld_comm_h
integer :: newprocs_comm_kgrid_except0,newrank_comm_kgrid_except0,newworld_comm_kgrid_except0
integer :: newprocs_comm_korbital_Vhxc,newrank_comm_korbital_Vhxc,newworld_comm_korbital_Vhxc
integer :: newprocs_bound(3),newrank_bound(3),new_world_bound(3)

integer :: num_pole_myrank
integer,allocatable :: icorr_polenum(:)
integer,allocatable :: icount_pole(:)
integer,allocatable :: icorr_xyz_pole(:,:,:)

integer,allocatable :: icoobox_bound(:,:,:)

! FFTE routine
integer :: ICOMMY,ICOMMZ
integer :: newprocs_ICOMMY,newprocs_ICOMMZ
integer :: newrank_ICOMMY,newrank_ICOMMZ
integer :: ICOMMY_copy,ICOMMZ_copy
integer :: newprocs_ICOMMY_copy,newprocs_ICOMMZ_copy
integer :: newrank_ICOMMY_copy,newrank_ICOMMZ_copy
integer :: iquot
integer :: i11,i12,i13,i14,i15
integer :: icheck_ascorder

CONTAINS

!=======================================================================
subroutine make_new_world
implicit none
integer :: i
integer :: i1,i2,i3,i4
integer :: ix,iy,iz
integer :: ixs,iys,izs
integer :: ibox
integer :: icolor,ikey

!only for identifying spin
!new_world for comm_spin
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i4=0,nproc_ob-1
      ibox=i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob
      if(myrank==ibox)then
        icolor=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob
        ikey=i4
      end if
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+i4*nproc_Mxin_mul
      if(myrank==ibox)then
        icolor=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))
        ikey=i4
      end if
    end do
  end do
  end do
  end do
end if

call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey,newworld_comm_spin,ierr)
call MPI_comm_size(newworld_comm_spin,newprocs_comm_spin,ierr)
call MPI_comm_rank(newworld_comm_spin,newrank_comm_spin,ierr)

!new_world for comm_grid
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i4=0,nproc_ob-1
      ibox=i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob
      if(ilsda==0)then
        if(myrank==ibox)then
          icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
          ikey=i4
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(myrank==ibox)then
            icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
            ikey=i4
          end if
        else
          if(myrank==ibox)then
            icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
            ikey=i4-nproc_ob_spin(1)
          end if
        end if
      end if
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+i4*nproc_Mxin_mul
      if(ilsda==0)then
        if(myrank==ibox)then
          icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
          ikey=i4
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(myrank==ibox)then
            icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
            ikey=i4
          end if
        else
          if(myrank==ibox)then
            icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
            ikey=i4-nproc_ob_spin(1)
          end if
        end if
      end if
    end do
  end do
  end do
  end do
end if
 
call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey,newworld_comm_grid,ierr)
call MPI_comm_size(newworld_comm_grid,newprocs_comm_grid,ierr)
call MPI_comm_rank(newworld_comm_grid,newrank_comm_grid,ierr)

!new_world for comm_orbital
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i4=0,nproc_ob-1
      ibox=i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob
      if(myrank==ibox)then
        icolor=i4
        ikey=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
      end if
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+i4*nproc_Mxin_mul
      if(myrank==ibox)then
        icolor=i4
        ikey=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
      end if
    end do
  end do
  end do
  end do
end if

call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey,newworld_comm_orbital,ierr)
call MPI_comm_size(newworld_comm_orbital,newprocs_comm_orbital,ierr)
call MPI_comm_rank(newworld_comm_orbital,newrank_comm_orbital,ierr)

!new_world for comm_mesh_s

nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

if(isequential==1)then
  do i4=0,nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    do i=0,nproc_Mxin_mul-1
      ibox=i1+i2*nproc_Mxin_s_dm(1)   &
             +i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
             +i4*nproc_Mxin_mul_s_dm  &
             +i*nproc/nproc_Mxin_mul
      if(myrank==ibox)then
        icolor=i4
        ikey=i1+i2*nproc_Mxin_s_dm(1)   &
               +i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
               +i*nproc_Mxin_mul_s_dm
      end if
    end do
  end do
  end do
  end do
  end do
else if(isequential==2)then
  do i4=0,nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    do i=0,nproc_Mxin_mul-1
      ibox=i+i1*nproc_Mxin_mul+i2*nproc_Mxin_mul*nproc_Mxin_s_dm(1)   &
            +i3*nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
            +i4*nproc_Mxin_mul*nproc_Mxin_mul_s_dm
      if(myrank==ibox)then
        icolor=i4
        ikey=i+i1*nproc_Mxin_mul+i2*nproc_Mxin_mul*nproc_Mxin_s_dm(1)   &
              +i3*nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
      end if
    end do
  end do
  end do
  end do
  end do
end if

call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey,newworld_comm_h,ierr)
call MPI_comm_size(newworld_comm_h,newprocs_comm_h,ierr)
call MPI_comm_rank(newworld_comm_h,newrank_comm_h,ierr)

if(isequential==1)then
  do iz=0,nproc_Mxin(3)-1
  do iy=0,nproc_Mxin(2)-1
  do ix=0,nproc_Mxin(1)-1
    do i4=0,nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
    do izs=0,nproc_Mxin_s_dm(3)-1
    do iys=0,nproc_Mxin_s_dm(2)-1
    do ixs=0,nproc_Mxin_s_dm(1)-1
      ibox=ixs+iys*nproc_Mxin_s_dm(1)   &
              +izs*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)   &
              +i4*nproc_Mxin_mul_s_dm    &
              +ix*nproc/nproc_Mxin_mul    &
              +iy*nproc/nproc_Mxin_mul*nproc_Mxin(1)   &
              +iz*nproc/nproc_Mxin_mul*nproc_Mxin(1)*nproc_Mxin(2)
      if(myrank==ibox)then
        imr(1)=ix
        imr(2)=iy
        imr(3)=iz
        imrs(1)=ixs
        imrs(2)=iys
        imrs(3)=izs
        igroup=i4
      end if
    end do 
    end do 
    end do 
  end do 
  end do 
  end do
  end do
else if(isequential==2)then
  do i4=0,nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
  do izs=0,nproc_Mxin_s_dm(3)-1
  do iys=0,nproc_Mxin_s_dm(2)-1
  do ixs=0,nproc_Mxin_s_dm(1)-1
    do iz=0,nproc_Mxin(3)-1
    do iy=0,nproc_Mxin(2)-1
    do ix=0,nproc_Mxin(1)-1
      ibox=ix+iy*nproc_Mxin(1)+iz*nproc_Mxin(1)*nproc_Mxin(2)  &
             +ixs*nproc_Mxin_mul  &
             +iys*nproc_Mxin_mul*nproc_Mxin_s_dm(1)  &
             +izs*nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
             +i4*nproc_Mxin_mul*nproc_Mxin_mul_s_dm
      if(myrank==ibox)then
        imr(1)=ix
        imr(2)=iy
        imr(3)=iz
        imrs(1)=ixs
        imrs(2)=iys
        imrs(3)=izs
        igroup=i4
      end if
    end do 
    end do 
    end do 
  end do 
  end do 
  end do
  end do
end if



if(isequential==1)then
  icolor=imrs(2)+imrs(3)*nproc_Mxin_s_dm(2)   &
                +igroup*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)   &
                +imr(2)*nproc/nproc_Mxin_mul/nproc_Mxin_s_dm(1)   &
                +imr(3)*nproc/nproc_Mxin_mul/nproc_Mxin_s_dm(1)*nproc_Mxin(2)
  ikey=imrs(1)+imr(1)*nproc_Mxin_s_dm(1)
else if(isequential==2)then
  icolor=imr(2)+imr(3)*nproc_Mxin(2)   &
               +imrs(2)*nproc_Mxin(2)*nproc_Mxin(3)   &
               +imrs(3)*nproc_Mxin(2)*nproc_Mxin(3)*nproc_Mxin_s_dm(2)  &
               +myrank/(nproc_Mxin_mul*nproc_Mxin_mul_s_dm)*nproc_Mxin(2)*nproc_Mxin(3)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
  ikey=imrs(1)+imr(1)*nproc_Mxin_s_dm(1)
end if

call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey,new_world_bound(1),ierr)
call MPI_comm_size(new_world_bound(1),newprocs_bound(1),ierr)
call MPI_comm_rank(new_world_bound(1),newrank_bound(1),ierr)

if(isequential==1)then
  icolor=imrs(1)+imrs(3)*nproc_Mxin_s_dm(1)   &
                +igroup*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(3)   &
                +imr(1)*nproc/nproc_Mxin_mul/nproc_Mxin_s_dm(2)   &
                +imr(3)*nproc/nproc_Mxin_mul/nproc_Mxin_s_dm(2)*nproc_Mxin(1)
  ikey=imrs(2)+imr(2)*nproc_Mxin_s_dm(2)
else if(isequential==2)then
  icolor=imr(1)+imr(3)*nproc_Mxin(1)   &
               +imrs(1)*nproc_Mxin(1)*nproc_Mxin(3)   &
               +imrs(3)*nproc_Mxin(1)*nproc_Mxin(3)*nproc_Mxin_s_dm(1)  &
               +myrank/(nproc_Mxin_mul*nproc_Mxin_mul_s_dm)*nproc_Mxin(1)*nproc_Mxin(3)*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(3)
  ikey=imrs(2)+imr(2)*nproc_Mxin_s_dm(2)
end if

call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey,new_world_bound(2),ierr)
call MPI_comm_size(new_world_bound(2),newprocs_bound(2),ierr)
call MPI_comm_rank(new_world_bound(2),newrank_bound(2),ierr)

if(isequential==1)then
  icolor=imrs(1)+imrs(2)*nproc_Mxin_s_dm(1)   &
                +igroup*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)   &
                +imr(1)*nproc/nproc_Mxin_mul/nproc_Mxin_s_dm(3)   &
                +imr(2)*nproc/nproc_Mxin_mul/nproc_Mxin_s_dm(3)*nproc_Mxin(1)
  ikey=imrs(3)+imr(3)*nproc_Mxin_s_dm(3)
else if(isequential==2)then
  icolor=imr(1)+imr(2)*nproc_Mxin(1)   &
               +imrs(1)*nproc_Mxin(1)*nproc_Mxin(2)   &
               +imrs(2)*nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin_s_dm(1)  &
               +myrank/(nproc_Mxin_mul*nproc_Mxin_mul_s_dm)*nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
  ikey=imrs(3)+imr(3)*nproc_Mxin_s_dm(3)
end if

call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey,new_world_bound(3),ierr)
call MPI_comm_size(new_world_bound(3),newprocs_bound(3),ierr)
call MPI_comm_rank(new_world_bound(3),newrank_bound(3),ierr)

if(isequential==1)then
  do i=0,nproc_Mxin_mul-1
    do i4=0,nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
    do i3=0,nproc_Mxin_s_dm(3)-1
    do i2=0,nproc_Mxin_s_dm(2)-1
    do i1=0,nproc_Mxin_s_dm(1)-1
      ibox=i1+i2*nproc_Mxin_s_dm(1)   &
             +i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)   &
             +i4*nproc_Mxin_mul_s_dm   &
             +i*nproc/nproc_Mxin_mul
      if(myrank==ibox)then
        icolor=i4+i*nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm
        ikey=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
      end if
    end do
    end do
    end do
    end do
  end do
else if(isequential==2)then
  do i4=0,nproc/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    do i=0,nproc_Mxin_mul-1
      ibox=i+i1*nproc_Mxin_mul+i2*nproc_Mxin_mul*nproc_Mxin_s_dm(1)   &
            +i3*nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
            +i4*nproc_Mxin_mul*nproc_Mxin_mul_s_dm
      if(myrank==ibox)then
        icolor=i+i4*nproc_Mxin_mul
        ikey=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
      end if
    end do
  end do
  end do
  end do
  end do
end if

call MPI_comm_split(MPI_COMM_WORLD,icolor,ikey,newworld_comm_korbital_Vhxc,ierr)
call MPI_comm_size(newworld_comm_h,newprocs_comm_korbital_Vhxc,ierr)
call MPI_comm_rank(newworld_comm_h,newrank_comm_korbital_Vhxc,ierr)

end subroutine make_new_world

!=====================================================================
subroutine make_corr_pole
implicit none

integer :: a,i
integer :: ix,iy,iz
integer :: ibox
integer :: j1,j2,j3
integer,allocatable :: ista_Mxin_pole(:,:)
integer,allocatable :: iend_Mxin_pole(:,:)
integer,allocatable :: inum_Mxin_pole(:,:)
integer,allocatable :: iflag_pole(:)
integer :: amin,amax
real(8) :: rmin,r
real(8),allocatable :: Rion2(:,:)
integer,allocatable :: nearatomnum(:,:,:)
integer,allocatable :: inv_icorr_polenum(:)

if(MEO==2)then

  if(iflag_ps==1)then
    amax=MI
    allocate(Rion2(3,MI))
    Rion2(:,:)=Rion(:,:)
  end if

  allocate(icount_pole(1:amax))
  allocate(nearatomnum(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
  icount_pole=0
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rmin=1.d6
    do a=1,amax
      r=sqrt( (gridcoo(ix,1)-Rion2(1,a))**2      &
            + (gridcoo(iy,2)-Rion2(2,a))**2      &
            + (gridcoo(iz,3)-Rion2(3,a))**2 )
      if ( r < rmin ) then
        rmin=r ; amin=a
      end if
    end do
    icount_pole(amin)=icount_pole(amin)+1
    nearatomnum(ix,iy,iz)=amin
  end do
  end do
  end do

  allocate(icorr_polenum(1:amax))
  allocate(inv_icorr_polenum(1:amax))
  icorr_polenum=0
  inv_icorr_polenum=0
  ibox=0
  do a=1,amax
    if(icount_pole(a)>=1)then
      ibox=ibox+1
      icorr_polenum(ibox)=a
      inv_icorr_polenum(a)=ibox
    end if
  end do
  num_pole_myrank=ibox

  allocate(icorr_xyz_pole(3,maxval(icount_pole(:)),num_pole_myrank))

  icount_pole(:)=0

  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    ibox=inv_icorr_polenum(nearatomnum(ix,iy,iz))
    icount_pole(ibox)=icount_pole(ibox)+1
    icorr_xyz_pole(1,icount_pole(ibox),ibox)=ix
    icorr_xyz_pole(2,icount_pole(ibox),ibox)=iy
    icorr_xyz_pole(3,icount_pole(ibox),ibox)=iz
  end do
  end do
  end do

  if(iflag_ps==1)then
    deallocate(Rion2)
  end if
  deallocate(nearatomnum)
  deallocate(inv_icorr_polenum)

else if(MEO==3)then

  allocate(ista_Mxin_pole(3,0:num_pole-1))
  allocate(iend_Mxin_pole(3,0:num_pole-1))
  allocate(inum_Mxin_pole(3,0:num_pole-1))
  allocate(iflag_pole(1:num_pole))

  do j3=0,num_pole_xyz(3)-1
  do j2=0,num_pole_xyz(2)-1
  do j1=0,num_pole_xyz(1)-1
    ibox = j1 + num_pole_xyz(1)*j2 + num_pole_xyz(1)*num_pole_xyz(2)*j3 
    ista_Mxin_pole(1,ibox)=j1*lg_num(1)/num_pole_xyz(1)+lg_sta(1)
    iend_Mxin_pole(1,ibox)=(j1+1)*lg_num(1)/num_pole_xyz(1)+lg_sta(1)-1
    ista_Mxin_pole(2,ibox)=j2*lg_num(2)/num_pole_xyz(2)+lg_sta(2)
    iend_Mxin_pole(2,ibox)=(j2+1)*lg_num(2)/num_pole_xyz(2)+lg_sta(2)-1
    ista_Mxin_pole(3,ibox)=j3*lg_num(3)/num_pole_xyz(3)+lg_sta(3)
    iend_Mxin_pole(3,ibox)=(j3+1)*lg_num(3)/num_pole_xyz(3)+lg_sta(3)-1
  end do
  end do
  end do

  iflag_pole=0

  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    do i=1,num_pole
      if(ista_Mxin_pole(3,i-1)<=iz.and.iend_Mxin_pole(3,i-1)>=iz.and.   &
         ista_Mxin_pole(2,i-1)<=iy.and.iend_Mxin_pole(2,i-1)>=iy.and.   &
         ista_Mxin_pole(1,i-1)<=ix.and.iend_Mxin_pole(1,i-1)>=ix)then
        iflag_pole(i)=1
      end if
    end do
  end do
  end do
  end do

  num_pole_myrank=0
  do i=1,num_pole
    if(iflag_pole(i)==1)then
      num_pole_myrank=num_pole_myrank+1
    end if
  end do

  allocate(icorr_polenum(1:num_pole_myrank))
  allocate(icount_pole(1:num_pole_myrank))

  ibox=1
  do i=1,num_pole
    if(iflag_pole(i)==1)then
      icorr_polenum(ibox)=i
      ibox=ibox+1
    end if
  end do

  icount_pole=0

  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    do i=1,num_pole_myrank
      if(ista_Mxin_pole(3,icorr_polenum(i)-1)<=iz.and.iend_Mxin_pole(3,icorr_polenum(i)-1)>=iz.and.   &
         ista_Mxin_pole(2,icorr_polenum(i)-1)<=iy.and.iend_Mxin_pole(2,icorr_polenum(i)-1)>=iy.and.   &
         ista_Mxin_pole(1,icorr_polenum(i)-1)<=ix.and.iend_Mxin_pole(1,icorr_polenum(i)-1)>=ix)then
        icount_pole(i)=icount_pole(i)+1
      end if
    end do
  end do
  end do
  end do
  
  allocate(icorr_xyz_pole(3,maxval(icount_pole(:)),num_pole_myrank))
 
  icount_pole=0

  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    do i=1,num_pole_myrank
      if(ista_Mxin_pole(3,icorr_polenum(i)-1)<=iz.and.iend_Mxin_pole(3,icorr_polenum(i)-1)>=iz.and.   &
         ista_Mxin_pole(2,icorr_polenum(i)-1)<=iy.and.iend_Mxin_pole(2,icorr_polenum(i)-1)>=iy.and.   &
         ista_Mxin_pole(1,icorr_polenum(i)-1)<=ix.and.iend_Mxin_pole(1,icorr_polenum(i)-1)>=ix)then
        icount_pole(i)=icount_pole(i)+1
        icorr_xyz_pole(1,icount_pole(i),i)=ix
        icorr_xyz_pole(2,icount_pole(i),i)=iy
        icorr_xyz_pole(3,icount_pole(i),i)=iz
      end if
    end do
  end do
  end do
  end do

  deallocate(ista_Mxin_pole,iend_Mxin_pole,inum_Mxin_pole)
  deallocate(iflag_pole)

end if

end subroutine make_corr_pole

!=====================================================================
subroutine make_icoobox_bound
implicit none
integer :: ix,iy,iz
integer :: ibox
integer :: icount

ibox=inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)*inum_Mxin_s(3,myrank)/minval(inum_Mxin_s(1:3,myrank))*2*Ndh

allocate( icoobox_bound(3,ibox,3) )

icount=0
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=lg_sta(1)-Ndh,lg_sta(1)-1
  icount=icount+1
  icoobox_bound(1,icount,1)=ix
  icoobox_bound(2,icount,1)=iy
  icoobox_bound(3,icount,1)=iz
end do
end do
end do
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=lg_end(1)+1,lg_end(1)+Ndh
  icount=icount+1
  icoobox_bound(1,icount,1)=ix
  icoobox_bound(2,icount,1)=iy
  icoobox_bound(3,icount,1)=iz
end do
end do
end do
icount=0
do iz=ng_sta(3),ng_end(3)
do iy=lg_sta(2)-Ndh,lg_sta(2)-1
do ix=ng_sta(1),ng_end(1)
  icount=icount+1
  icoobox_bound(1,icount,2)=ix
  icoobox_bound(2,icount,2)=iy
  icoobox_bound(3,icount,2)=iz
end do
end do
end do
do iz=ng_sta(3),ng_end(3)
do iy=lg_end(2)+1,lg_end(2)+Ndh
do ix=ng_sta(1),ng_end(1)
  icount=icount+1
  icoobox_bound(1,icount,2)=ix
  icoobox_bound(2,icount,2)=iy
  icoobox_bound(3,icount,2)=iz
end do
end do
end do
icount=0
do iz=lg_sta(3)-Ndh,lg_sta(3)-1
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  icount=icount+1
  icoobox_bound(1,icount,3)=ix
  icoobox_bound(2,icount,3)=iy
  icoobox_bound(3,icount,3)=iz
end do
end do
end do
do iz=lg_end(3)+1,lg_end(3)+Ndh
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  icount=icount+1
  icoobox_bound(1,icount,3)=ix
  icoobox_bound(2,icount,3)=iy
  icoobox_bound(3,icount,3)=iz
end do
end do
end do

return
end subroutine make_icoobox_bound

!=====================================================================
!======================================================================
subroutine mpi_allgatherv_vlocal
!$ use omp_lib

implicit none
integer :: i
integer :: i1,i2,i3
integer :: ix,iy,iz
integer :: ibox,ibox2,ibox3
real(8),allocatable :: matbox11(:),matbox12(:)
integer :: iscnt
integer,allocatable :: ircnt(:)
integer,allocatable :: idisp(:)
integer :: is,is_sta,is_end

elp3(1001)=MPI_Wtime()

elp3(1002)=MPI_Wtime()
elp3(1052)=elp3(1052)+elp3(1002)-elp3(1001)

allocate(ircnt(0:nproc_Mxin_mul_s_dm-1))
allocate(idisp(0:nproc_Mxin_mul_s_dm-1))

allocate (matbox11(0:(inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)*inum_Mxin_s(3,myrank))-1))
allocate (matbox12(0:(mg_num(1)*mg_num(2)*mg_num(3))-1))

iscnt=inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)*inum_Mxin_s(3,myrank)
if(isequential==1)then
  do i=0,nproc_Mxin_mul_s_dm-1
    ibox=(myrank/nproc_Mxin_mul_s_dm)*nproc_Mxin_mul_s_dm+i
    ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
  end do
else if(isequential==2)then
  do i=0,nproc_Mxin_mul_s_dm-1
    ibox=mod(myrank,nproc_Mxin_mul)+i*nproc_Mxin_mul
    ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
  end do
end if

idisp(0)=0
do i=1,nproc_Mxin_mul_s_dm-1
  idisp(i)=idisp(i-1)+ircnt(i-1)
end do

if(ilsda==0)then
  is_sta=1
  is_end=1
else
  is_sta=1
  is_end=2
end if

do is=is_sta,is_end

  if(iSCFRT==1)then
!$OMP parallel do private(ibox3,ix,iy,iz)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,myrank)   &
                    +(iz-ng_sta(3))*inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)
      matbox11(ibox3) = Vpsl(ix,iy,iz)+Vh(ix,iy,iz)
    end do
    end do
    end do
  else if(iSCFRT==2)then
    if(mod(itt,2)==1)then
!$OMP parallel do private(ibox3,ix,iy,iz)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,myrank)   &
                      +(iz-ng_sta(3))*inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)
        matbox11(ibox3) = Vpsl(ix,iy,iz)+Vh_stock2(ix,iy,iz)
      end do
      end do
      end do
    else
!$OMP parallel do private(ibox3,ix,iy,iz)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,myrank)   &
                      +(iz-ng_sta(3))*inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)
        matbox11(ibox3) = Vpsl(ix,iy,iz)+Vh_stock1(ix,iy,iz)
      end do
      end do
      end do
    end if
  end if
  if(ilsda==0)then
!$OMP parallel do private(ibox3,ix,iy,iz)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,myrank)   &
                    +(iz-ng_sta(3))*inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)
      matbox11(ibox3) = matbox11(ibox3)+Vxc(ix,iy,iz)
    end do
    end do
    end do
  else
!$OMP parallel do private(ibox3,ix,iy,iz)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,myrank)   &
                    +(iz-ng_sta(3))*inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)
      matbox11(ibox3) = matbox11(ibox3)+Vxc_s(ix,iy,iz,is)
    end do
    end do
    end do
  end if

  elp3(761)=MPI_Wtime()
  call MPI_Allgatherv(matbox11,iscnt,      MPI_DOUBLE_PRECISION,      &
                      matbox12,ircnt,idisp,MPI_DOUBLE_PRECISION,      &
                      newworld_comm_grid,ierr)
  elp3(762)=MPI_Wtime()
  elp3(781)=elp3(781)+elp3(762)-elp3(761) 


  if(isequential==1)then
!$OMP parallel do private(i1,i2,i3,ibox,ibox2) collapse(3)
    do i3=0,nproc_Mxin_s_dm(3)-1
    do i2=0,nproc_Mxin_s_dm(2)-1
    do i1=0,nproc_Mxin_s_dm(1)-1
      ibox=(myrank/nproc_Mxin_mul_s_dm)*nproc_Mxin_mul_s_dm    &
            +(i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2))
      ibox2=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)

      call copyVlocal(matbox12(idisp(ibox2):  &
                      (idisp(ibox2)+inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox))),  &
                      ibox,ibox2,is)

    end do
    end do
    end do
  else if(isequential==2)then
!$OMP parallel do private(i1,i2,i3,ibox,ibox2) collapse(3)
    do i3=0,nproc_Mxin_s_dm(3)-1
    do i2=0,nproc_Mxin_s_dm(2)-1
    do i1=0,nproc_Mxin_s_dm(1)-1
      ibox=mod(myrank,nproc_Mxin_mul)    &
          +(i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2))*nproc_Mxin_mul
      ibox2=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)

      call copyVlocal(matbox12(idisp(ibox2):  &
                      (idisp(ibox2)+inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox))),  &
                      ibox,ibox2,is)

    end do
    end do
    end do
  end if

end do

deallocate (ircnt,idisp)
deallocate (matbox11)
deallocate (matbox12)

elp3(1003)=MPI_Wtime()
elp3(1053)=elp3(1053)+elp3(1003)-elp3(1002)
elp3(1054)=elp3(1054)+elp3(1003)-elp3(1001)

end subroutine mpi_allgatherv_vlocal

!======================================================================
subroutine mpibcast_mesh_s_kxc(Vbox)
!$ use omp_lib

implicit none
integer :: i
integer :: i1,i2,i3
integer :: ix,iy,iz
integer :: ibox,ibox2,ibox3
real(8),allocatable :: matbox1(:),matbox2(:)
real(8) :: Vbox(mg_sta(1):mg_end(1),  &
                mg_sta(2):mg_end(2),  &
                mg_sta(3):mg_end(3))
integer :: iscnt
integer,allocatable :: ircnt(:)
integer,allocatable :: idisp(:)

if(nproc_ob/=1.or.nproc_Mxin_mul/=1)then

  allocate(ircnt(0:nproc_Mxin_mul_s_dm-1))
  allocate(idisp(0:nproc_Mxin_mul_s_dm-1))

  iscnt=inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)*inum_Mxin_s(3,myrank)
  do i=0,nproc_Mxin_mul_s_dm-1
    ibox=mod(myrank,nproc_Mxin_mul)+i*nproc_Mxin_mul
    ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
  end do

  idisp(0)=0
  do i=1,nproc_Mxin_mul_s_dm-1
    idisp(i)=idisp(i-1)+ircnt(i-1)
  end do

  allocate (matbox1(0:(inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)*inum_Mxin_s(3,myrank))-1))
  allocate (matbox2(0:(mg_num(1)*mg_num(2)*mg_num(3))-1))
  
!$OMP parallel do private(ibox3,ix,iy,iz)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,myrank)   &
                  +(iz-ng_sta(3))*inum_Mxin_s(1,myrank)*inum_Mxin_s(2,myrank)
    matbox1(ibox3) = Vbox(ix,iy,iz)
  end do
  end do
  end do

  call MPI_Allgatherv(matbox1,iscnt,      MPI_DOUBLE_PRECISION,      &
                      matbox2,ircnt,idisp,MPI_DOUBLE_PRECISION,      &
                      newworld_comm_h,ierr)
   
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    ibox=mod(myrank,nproc_Mxin_mul)    &
          +(i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2))*nproc_Mxin_mul
    ibox2=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)

    if(newrank_comm_h/=ibox2)then
!$OMP parallel do private(ibox3,ix,iy,iz)
      do iz=ista_Mxin_s(3,ibox),iend_Mxin_s(3,ibox)
      do iy=ista_Mxin_s(2,ibox),iend_Mxin_s(2,ibox)
      do ix=ista_Mxin_s(1,ibox),iend_Mxin_s(1,ibox)
        ibox3=idisp(ibox2)+ix-ista_Mxin_s(1,ibox)+(iy-ista_Mxin_s(2,ibox))*inum_Mxin_s(1,ibox)   &
                      +(iz-ista_Mxin_s(3,ibox))*inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)
        Vbox(ix,iy,iz) = matbox2(ibox3) 
      end do
      end do
      end do
    end if
  end do
  end do
  end do

  deallocate (ircnt,idisp)
  deallocate (matbox1)
  deallocate (matbox2)

end if

end subroutine mpibcast_mesh_s_kxc

END MODULE new_world_sub

