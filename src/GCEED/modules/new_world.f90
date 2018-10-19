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

integer :: num_pole_myrank
integer,allocatable :: icorr_polenum(:)
integer,allocatable :: icount_pole(:)
integer,allocatable :: icorr_xyz_pole(:,:,:)

integer,allocatable :: icoobox_bound(:,:,:)

! FFTE routine
integer :: iquot
integer :: i11,i12,i13,i14,i15
integer :: icheck_ascorder

CONTAINS

!=======================================================================
subroutine make_new_world
use salmon_parallel
use salmon_communication, only: comm_create_group, comm_get_groupinfo, &
                                comm_summation
use misc_routines, only: get_wtime
implicit none
integer :: ii,jj
integer :: i1,i2,i3,i4,i5
integer :: ix,iy,iz
integer :: ixs,iys,izs
integer :: ibox
integer :: icolor,ikey

integer :: iarray_ascorder(0:nproc_size_global-1)
integer :: iarray_ascorder2(0:nproc_size_global-1)
integer :: icheck_ascorder_tmp
integer :: icheck_ascorder_tmp2
integer,allocatable :: iwork(:),iwork2(:)

integer :: LNPU2(3)

!new_world for comm_kgrid
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob*nproc_k
      if(nproc_id_global==ibox)then
        icolor=i5*nproc_ob+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob*nproc_k
        ikey=i4
      end if
    end do
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+(i5*nproc_ob+i4)*nproc_Mxin_mul
      if(nproc_id_global==ibox)then
        icolor=i5+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_k
        ikey=i4
      end if
    end do
    end do
  end do
  end do
  end do
end if

nproc_group_kgrid = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_kgrid, nproc_id_kgrid, nproc_size_kgrid)

!only for identifying spin
!new_world for comm_spin
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob*nproc_k
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i5+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_k
          ikey=i4
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=0+2*(i5+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_k)
            ikey=i4
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=1+2*(i5+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_k)
            ikey=i4-nproc_ob_spin(1)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+(i5*nproc_ob+i4)*nproc_Mxin_mul
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i5+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_k
          ikey=i4
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=0+2*(i5+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_k)
            ikey=i4
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=1+2*(i5+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_k)
            ikey=i4-nproc_ob_spin(1)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
end if
 
nproc_group_spin = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_spin, nproc_id_spin, nproc_size_spin)

!new_world for comm_korbital
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob*nproc_k
      if(nproc_id_global==ibox)then
        icolor=i5*nproc_ob+i4
        ikey=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
      end if
    end do
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+(i5*nproc_ob+i4)*nproc_Mxin_mul
      if(nproc_id_global==ibox)then
        icolor=i5*nproc_ob+i4
        ikey=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
      end if
    end do
    end do
  end do
  end do
  end do
end if

nproc_group_korbital = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_korbital, nproc_id_korbital, nproc_size_korbital)

!new_world for comm_rho
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob*nproc_k
      if(nproc_id_global==ibox)then
        icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
        ikey=i5*nproc_ob+i4
      end if
    end do
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+(i5*nproc_ob+i4)*nproc_Mxin_mul
      if(nproc_id_global==ibox)then
        icolor=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))
        ikey=i5*nproc_ob+i4
      end if
    end do
    end do
  end do
  end do
  end do
end if
 
nproc_group_rho = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_rho, nproc_id_rho, nproc_size_rho)

!new_world for comm_k
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob*nproc_k
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i5
          ikey=i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=2*i5+0
            ikey=i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob_spin(1)
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=2*i5+1
            ikey=i4-nproc_ob_spin(1)+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob_spin(2)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+(i5*nproc_ob+i4)*nproc_Mxin_mul
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i5
          ikey=i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=2*i5+0
            ikey=i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob_spin(1)
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=2*i5+1
            ikey=i4-nproc_ob_spin(1)+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob_spin(2)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
end if
 
nproc_group_k = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_k, nproc_id_k, nproc_size_k)

!new_world for comm_grid
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob*nproc_k
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
          ikey=i5*nproc_ob+i4
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
            ikey=i5*nproc_ob_spin(1)+i4
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
            ikey=i5*nproc_ob_spin(2)+i4-nproc_ob_spin(1)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+(i5*nproc_ob+i4)*nproc_Mxin_mul
      if(ilsda==0)then
        if(nproc_id_global==ibox)then
          icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
          ikey=i5*nproc_ob+i4
        end if
      else
        if(i4<nproc_ob_spin(1))then
          if(nproc_id_global==ibox)then
            icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
            ikey=i5*nproc_ob_spin(1)+i4
          end if
        else
          if(nproc_id_global==ibox)then
            icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)
            ikey=i5*nproc_ob_spin(2)+i4-nproc_ob_spin(1)
          end if
        end if
      end if
    end do
    end do
  end do
  end do
  end do
end if
 
nproc_group_grid = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_grid, nproc_id_grid, nproc_size_grid)

!new_world for comm_orbitalgrid
if(isequential==1)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=i5*nproc_ob+i4+(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))*nproc_ob*nproc_k
      if(nproc_id_global==ibox)then
        icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)+i4*nproc_Mxin_mul
        ikey=i5
      end if
    end do
    end do
  end do
  end do
  end do
else if(isequential==2)then
  do i3=0,nproc_Mxin(3)-1
  do i2=0,nproc_Mxin(2)-1
  do i1=0,nproc_Mxin(1)-1
    do i5=0,nproc_k-1
    do i4=0,nproc_ob-1
      ibox=(i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2))+(i5*nproc_ob+i4)*nproc_Mxin_mul
      if(nproc_id_global==ibox)then
        icolor=i1+i2*nproc_Mxin(1)+i3*nproc_Mxin(1)*nproc_Mxin(2)+i4*nproc_Mxin_mul
        ikey=i5
      end if
    end do
    end do
  end do
  end do
  end do
end if
 
nproc_group_orbitalgrid = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_orbitalgrid, nproc_id_orbitalgrid, nproc_size_orbitalgrid)

!new_world for comm_mesh_s

nproc_Mxin_mul_s_dm=nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)

if(isequential==1)then
  do i4=0,nproc_size_global/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    do ii=0,nproc_Mxin_mul-1
      ibox=i1+i2*nproc_Mxin_s_dm(1)   &
             +i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
             +i4*nproc_Mxin_mul_s_dm  &
             +ii*nproc_size_global/nproc_Mxin_mul
      if(nproc_id_global==ibox)then
        icolor=i4
        ikey=i1+i2*nproc_Mxin_s_dm(1)   &
               +i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
               +ii*nproc_Mxin_mul_s_dm
      end if
    end do
  end do
  end do
  end do
  end do
else if(isequential==2)then
  do i4=0,nproc_size_global/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    do ii=0,nproc_Mxin_mul-1
      ibox=ii+i1*nproc_Mxin_mul+i2*nproc_Mxin_mul*nproc_Mxin_s_dm(1)   &
            +i3*nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
            +i4*nproc_Mxin_mul*nproc_Mxin_mul_s_dm
      if(nproc_id_global==ibox)then
        icolor=i4
        ikey=ii+i1*nproc_Mxin_mul+i2*nproc_Mxin_mul*nproc_Mxin_s_dm(1)   &
              +i3*nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
      end if
    end do
  end do
  end do
  end do
  end do
end if

nproc_group_h = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_h, nproc_id_h, nproc_size_h)

if(isequential==1)then
  do iz=0,nproc_Mxin(3)-1
  do iy=0,nproc_Mxin(2)-1
  do ix=0,nproc_Mxin(1)-1
    do i4=0,nproc_size_global/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
    do izs=0,nproc_Mxin_s_dm(3)-1
    do iys=0,nproc_Mxin_s_dm(2)-1
    do ixs=0,nproc_Mxin_s_dm(1)-1
      ibox=ixs+iys*nproc_Mxin_s_dm(1)   &
              +izs*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)   &
              +i4*nproc_Mxin_mul_s_dm    &
              +ix*nproc_size_global/nproc_Mxin_mul    &
              +iy*nproc_size_global/nproc_Mxin_mul*nproc_Mxin(1)   &
              +iz*nproc_size_global/nproc_Mxin_mul*nproc_Mxin(1)*nproc_Mxin(2)
      if(nproc_id_global==ibox)then
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
  do i4=0,nproc_size_global/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
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
      if(nproc_id_global==ibox)then
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
                +imr(2)*nproc_size_global/nproc_Mxin_mul/nproc_Mxin_s_dm(1)   &
                +imr(3)*nproc_size_global/nproc_Mxin_mul/nproc_Mxin_s_dm(1)*nproc_Mxin(2)
  ikey=imrs(1)+imr(1)*nproc_Mxin_s_dm(1)
else if(isequential==2)then
  icolor=imr(2)+imr(3)*nproc_Mxin(2)   &
               +imrs(2)*nproc_Mxin(2)*nproc_Mxin(3)   &
               +imrs(3)*nproc_Mxin(2)*nproc_Mxin(3)*nproc_Mxin_s_dm(2)  &
               +nproc_id_global/(nproc_Mxin_mul*nproc_Mxin_mul_s_dm)*nproc_Mxin(2)  &
               *nproc_Mxin(3)*nproc_Mxin_s_dm(2)*nproc_Mxin_s_dm(3)
!  ikey=imr(1)+imrs(1)*nproc_Mxin(1)
  ikey=imrs(1)+imr(1)*nproc_Mxin_s_dm(1)
end if

nproc_group_bound(1) = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_bound(1), nproc_id_bound(1), nproc_size_bound(1))

if(isequential==1)then
  icolor=imrs(1)+imrs(3)*nproc_Mxin_s_dm(1)   &
                +igroup*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(3)   &
                +imr(1)*nproc_size_global/nproc_Mxin_mul/nproc_Mxin_s_dm(2)   &
                +imr(3)*nproc_size_global/nproc_Mxin_mul/nproc_Mxin_s_dm(2)*nproc_Mxin(1)
  ikey=imrs(2)+imr(2)*nproc_Mxin_s_dm(2)
else if(isequential==2)then
  icolor=imr(1)+imr(3)*nproc_Mxin(1)   &
               +imrs(1)*nproc_Mxin(1)*nproc_Mxin(3)   &
               +imrs(3)*nproc_Mxin(1)*nproc_Mxin(3)*nproc_Mxin_s_dm(1)  &
               +nproc_id_global/(nproc_Mxin_mul*nproc_Mxin_mul_s_dm)*nproc_Mxin(1)*nproc_Mxin(3)  &
               *nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(3)
!  ikey=imr(2)+imrs(2)*nproc_Mxin(2)
  ikey=imrs(2)+imr(2)*nproc_Mxin_s_dm(2)
end if

nproc_group_bound(2) = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_bound(2), nproc_id_bound(2), nproc_size_bound(2))

if(isequential==1)then
  icolor=imrs(1)+imrs(2)*nproc_Mxin_s_dm(1)   &
                +igroup*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)   &
                +imr(1)*nproc_size_global/nproc_Mxin_mul/nproc_Mxin_s_dm(3)   &
                +imr(2)*nproc_size_global/nproc_Mxin_mul/nproc_Mxin_s_dm(3)*nproc_Mxin(1)
  ikey=imrs(3)+imr(3)*nproc_Mxin_s_dm(3)
else if(isequential==2)then
  icolor=imr(1)+imr(2)*nproc_Mxin(1)   &
               +imrs(1)*nproc_Mxin(1)*nproc_Mxin(2)   &
               +imrs(2)*nproc_Mxin(1)*nproc_Mxin(2)*nproc_Mxin_s_dm(1)  &
               +nproc_id_global/(nproc_Mxin_mul*nproc_Mxin_mul_s_dm)*nproc_Mxin(1)*nproc_Mxin(2)  &
               *nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
!  ikey=imr(3)+imrs(3)*nproc_Mxin(3)
  ikey=imrs(3)+imr(3)*nproc_Mxin_s_dm(3)
end if

nproc_group_bound(3) = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_bound(3), nproc_id_bound(3), nproc_size_bound(3))

if(isequential==1)then
  do ii=0,nproc_Mxin_mul-1
    do i4=0,nproc_size_global/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
    do i3=0,nproc_Mxin_s_dm(3)-1
    do i2=0,nproc_Mxin_s_dm(2)-1
    do i1=0,nproc_Mxin_s_dm(1)-1
      ibox=i1+i2*nproc_Mxin_s_dm(1)   &
             +i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)   &
             +i4*nproc_Mxin_mul_s_dm   &
             +ii*nproc_size_global/nproc_Mxin_mul
      if(nproc_id_global==ibox)then
        icolor=i4+ii*nproc_size_global/nproc_Mxin_mul/nproc_Mxin_mul_s_dm
        ikey=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
      end if
    end do
    end do
    end do
    end do
  end do
else if(isequential==2)then
  do i4=0,nproc_size_global/nproc_Mxin_mul/nproc_Mxin_mul_s_dm-1
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    do ii=0,nproc_Mxin_mul-1
      ibox=ii+i1*nproc_Mxin_mul+i2*nproc_Mxin_mul*nproc_Mxin_s_dm(1)   &
            +i3*nproc_Mxin_mul*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)  &
            +i4*nproc_Mxin_mul*nproc_Mxin_mul_s_dm
      if(nproc_id_global==ibox)then
        icolor=ii+i4*nproc_Mxin_mul
        ikey=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)
      end if
    end do
  end do
  end do
  end do
  end do
end if

nproc_group_korbital_vhxc = comm_create_group(nproc_group_global, icolor, ikey)
call comm_get_groupinfo(nproc_group_korbital_vhxc, nproc_id_korbital_vhxc, nproc_size_korbital_vhxc)

!call FACTOR(nproc,LNPU)
!NPUZ=(2**(LNPU(1)/2))*(3**(LNPU(2)/2))*(5**(LNPU(3)/2))
!NPUY=nproc/NPUZ

!if(iflag_hartree==4)then

  icheck_ascorder=1
  
  iarray_ascorder=0
  iarray_ascorder(nproc_id_bound(2))=nproc_id_global
  allocate(iwork(0:nproc_size_bound(2)-1),iwork2(0:nproc_size_bound(2)-1))
  iwork(0:nproc_size_bound(2)-1)=iarray_ascorder(0:nproc_size_bound(2)-1)
  call comm_summation(iwork,iwork2,nproc_size_bound(2),nproc_group_bound(2))
  iarray_ascorder2(0:nproc_size_bound(2)-1)=iwork2(0:nproc_size_bound(2)-1)
  deallocate(iwork,iwork2)

  icheck_ascorder_tmp=1
  do jj=0,nproc_size_bound(2)-2
    if(iarray_ascorder2(jj)>iarray_ascorder2(jj+1)) icheck_ascorder_tmp=0
  end do
  
  call comm_summation(icheck_ascorder_tmp,icheck_ascorder_tmp2,nproc_group_global)
 
  if(icheck_ascorder_tmp2/=nproc_size_global) icheck_ascorder=0

  iarray_ascorder=0
  iarray_ascorder(nproc_id_bound(3))=nproc_id_global
  allocate(iwork(0:nproc_size_bound(3)-1),iwork2(0:nproc_size_bound(3)-1))
  iwork(0:nproc_size_bound(3)-1)=iarray_ascorder(0:nproc_size_bound(3)-1)
  call comm_summation(iwork,iwork2,nproc_size_bound(3),nproc_group_bound(3))
  iarray_ascorder2(0:nproc_size_bound(3)-1)=iwork2(0:nproc_size_bound(3)-1)
  deallocate(iwork,iwork2)

  icheck_ascorder_tmp=1
  do jj=0,nproc_size_bound(3)-2
    if(iarray_ascorder2(jj)>iarray_ascorder2(jj+1)) icheck_ascorder_tmp=0
  end do
  
  call comm_summation(icheck_ascorder_tmp,icheck_ascorder_tmp2,nproc_group_global)
 
  if(icheck_ascorder_tmp2/=nproc_size_global) icheck_ascorder=0

  if(iperiodic==3.and.iflag_hartree==4)then
    if(nproc_id_global==0)then
      if(icheck_ascorder==0)then
        write(*,*) "Ranks are NOT in ascending order. Allreduce is done in FFTE routine."
      else if(icheck_ascorder==1)then
        write(*,*) "Ranks are in ascending order. Allreduce is skipped in FFTE routine."
      end if
    end if
  end if

! communicators for FFTE routine
  if(icheck_ascorder==1)then
    NPUW=nproc_Mxin_s_dm(1)*nproc_Mxin(1)
    NPUY=nproc_Mxin_s_dm(2)*nproc_Mxin(2)
    NPUZ=nproc_Mxin_s_dm(3)*nproc_Mxin(3)
 
    icolor=nproc_id_bound(3)+nproc_id_bound(1)*NPUZ 
    ikey=nproc_id_bound(2)
    nproc_group_icommy = comm_create_group(nproc_group_global, icolor, ikey)
    call comm_get_groupinfo(nproc_group_icommy, nproc_id_icommy, nproc_size_icommy)

    icolor=nproc_id_bound(2)+nproc_id_bound(1)*NPUY 
    ikey=nproc_id_bound(3)
    nproc_group_icommz = comm_create_group(nproc_group_global, icolor, ikey)
    call comm_get_groupinfo(nproc_group_icommz, nproc_id_icommz, nproc_size_icommz)

    icolor=nproc_id_bound(2)+nproc_id_bound(3)*NPUY 
    ikey=nproc_id_bound(1)
    nproc_group_icommw = comm_create_group(nproc_group_global, icolor, ikey)
    call comm_get_groupinfo(nproc_group_icommw, nproc_id_icommw, nproc_size_icommw)

  else
    call factor(nproc_size_global,LNPU2)
    NPUZ=(2**(LNPU2(1)/2))*(3**(LNPU2(2)/2))*(5**(LNPU2(3)/2))
    NPUY=nproc_size_global/NPUZ
    NPUW=1
  
    icolor=nproc_id_global/NPUY
    ikey=0
    nproc_group_icommy = comm_create_group(nproc_group_global, icolor, ikey)
    call comm_get_groupinfo(nproc_group_icommy, nproc_id_icommy, nproc_size_icommy)

    icolor=mod(nproc_id_global,NPUY)
    ikey=0
    nproc_group_icommz = comm_create_group(nproc_group_global, icolor, ikey)
    call comm_get_groupinfo(nproc_group_icommz, nproc_id_icommz, nproc_size_icommz)
    
  end if

  iquot=nproc_id_global/(NPUY*NPUZ)
  
  i11=mod(nproc_id_global,nproc_Mxin(2)*nproc_Mxin(3))
  i12=i11/nproc_Mxin(2)
  i13=i12*nproc_Mxin(3)
  i14=nproc_id_global/(NPUY*nproc_Mxin(3))
  icolor=i13+i14+iquot*NPUZ
  
  i11=mod(nproc_id_global,nproc_Mxin(2))
  i12=nproc_id_global/(nproc_Mxin(2)*nproc_Mxin(3))
  ikey=i11*NPUY/nproc_Mxin(2)+mod(i12,NPUY/nproc_Mxin(2))
  
  nproc_group_icommy_copy = comm_create_group(nproc_group_global, icolor, ikey)
  call comm_get_groupinfo(nproc_group_icommy_copy, nproc_id_icommy_copy, nproc_size_icommy_copy)

!end if



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
use salmon_parallel
implicit none
integer :: ix,iy,iz
integer :: ibox
integer :: icount

ibox=inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)  &
     *inum_Mxin_s(3,nproc_id_global)/minval(inum_Mxin_s(1:3,nproc_id_global))*2*Ndh

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
subroutine allgatherv_vlocal
use salmon_parallel, only: nproc_id_global, nproc_group_grid
use salmon_communication, only: comm_allgatherv
use misc_routines, only: get_wtime

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

elp3(1001)=get_wtime()

elp3(1002)=get_wtime()
elp3(1052)=elp3(1052)+elp3(1002)-elp3(1001)

allocate(ircnt(0:nproc_Mxin_mul_s_dm-1))
allocate(idisp(0:nproc_Mxin_mul_s_dm-1))

allocate (matbox11(0:(inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)*inum_Mxin_s(3,nproc_id_global))-1))
allocate (matbox12(0:(mg_num(1)*mg_num(2)*mg_num(3))-1))

iscnt=inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)*inum_Mxin_s(3,nproc_id_global)
if(isequential==1)then
  do i=0,nproc_Mxin_mul_s_dm-1
    ibox=(nproc_id_global/nproc_Mxin_mul_s_dm)*nproc_Mxin_mul_s_dm+i
    ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
  end do
else if(isequential==2)then
  do i=0,nproc_Mxin_mul_s_dm-1
    ibox=mod(nproc_id_global,nproc_Mxin_mul)+i*nproc_Mxin_mul
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
      ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,nproc_id_global)   &
                    +(iz-ng_sta(3))*inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)
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
        ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,nproc_id_global)   &
                      +(iz-ng_sta(3))*inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)
        matbox11(ibox3) = Vpsl(ix,iy,iz)+Vh_stock2(ix,iy,iz)
      end do
      end do
      end do
    else
!$OMP parallel do private(ibox3,ix,iy,iz)
      do iz=ng_sta(3),ng_end(3)
      do iy=ng_sta(2),ng_end(2)
      do ix=ng_sta(1),ng_end(1)
        ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,nproc_id_global)   &
                      +(iz-ng_sta(3))*inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)
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
      ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,nproc_id_global)   &
                    +(iz-ng_sta(3))*inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)
      matbox11(ibox3) = matbox11(ibox3)+Vxc(ix,iy,iz)
    end do
    end do
    end do
  else
!$OMP parallel do private(ibox3,ix,iy,iz)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,nproc_id_global)   &
                    +(iz-ng_sta(3))*inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)
      matbox11(ibox3) = matbox11(ibox3)+Vxc_s(ix,iy,iz,is)
    end do
    end do
    end do
  end if

  elp3(761)=get_wtime()
  call comm_allgatherv(matbox11,matbox12,ircnt,idisp,nproc_group_grid)
  elp3(762)=get_wtime()
  elp3(781)=elp3(781)+elp3(762)-elp3(761) 


  if(isequential==1)then
!$OMP parallel do private(i1,i2,i3,ibox,ibox2) collapse(3)
    do i3=0,nproc_Mxin_s_dm(3)-1
    do i2=0,nproc_Mxin_s_dm(2)-1
    do i1=0,nproc_Mxin_s_dm(1)-1
      ibox=(nproc_id_global/nproc_Mxin_mul_s_dm)*nproc_Mxin_mul_s_dm    &
            +(i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2))
      ibox2=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)

      call copyVlocal(matbox12(idisp(ibox2):  &
                      (idisp(ibox2)+inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox))),  &
                      ibox,is)

    end do
    end do
    end do
  else if(isequential==2)then
!$OMP parallel do private(i1,i2,i3,ibox,ibox2) collapse(3)
    do i3=0,nproc_Mxin_s_dm(3)-1
    do i2=0,nproc_Mxin_s_dm(2)-1
    do i1=0,nproc_Mxin_s_dm(1)-1
      ibox=mod(nproc_id_global,nproc_Mxin_mul)    &
          +(i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2))*nproc_Mxin_mul
      ibox2=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)

      call copyVlocal(matbox12(idisp(ibox2):  &
                      (idisp(ibox2)+inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox))),  &
                      ibox,is)

    end do
    end do
    end do
  end if

end do

deallocate (ircnt,idisp)
deallocate (matbox11)
deallocate (matbox12)

elp3(1003)=get_wtime()
elp3(1053)=elp3(1053)+elp3(1003)-elp3(1002)
elp3(1054)=elp3(1054)+elp3(1003)-elp3(1001)

end subroutine allgatherv_vlocal

!======================================================================
subroutine mpibcast_mesh_s_kxc(Vbox)
use salmon_parallel, only: nproc_id_global, nproc_group_h, nproc_id_h
use salmon_communication, only: comm_allgatherv

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

  iscnt=inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)*inum_Mxin_s(3,nproc_id_global)
  do i=0,nproc_Mxin_mul_s_dm-1
    ibox=mod(nproc_id_global,nproc_Mxin_mul)+i*nproc_Mxin_mul
    ircnt(i)=inum_Mxin_s(1,ibox)*inum_Mxin_s(2,ibox)*inum_Mxin_s(3,ibox)
  end do

  idisp(0)=0
  do i=1,nproc_Mxin_mul_s_dm-1
    idisp(i)=idisp(i-1)+ircnt(i-1)
  end do

  allocate (matbox1(0:(inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)*inum_Mxin_s(3,nproc_id_global))-1))
  allocate (matbox2(0:(mg_num(1)*mg_num(2)*mg_num(3))-1))
  
!$OMP parallel do private(ibox3,ix,iy,iz)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    ibox3=ix-ng_sta(1)+(iy-ng_sta(2))*inum_Mxin_s(1,nproc_id_global)   &
                  +(iz-ng_sta(3))*inum_Mxin_s(1,nproc_id_global)*inum_Mxin_s(2,nproc_id_global)
    matbox1(ibox3) = Vbox(ix,iy,iz)
  end do
  end do
  end do

  call comm_allgatherv(matbox1,matbox2,ircnt,idisp,nproc_group_h)
   
  do i3=0,nproc_Mxin_s_dm(3)-1
  do i2=0,nproc_Mxin_s_dm(2)-1
  do i1=0,nproc_Mxin_s_dm(1)-1
    ibox=mod(nproc_id_global,nproc_Mxin_mul)    &
          +(i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2))*nproc_Mxin_mul
    ibox2=i1+i2*nproc_Mxin_s_dm(1)+i3*nproc_Mxin_s_dm(1)*nproc_Mxin_s_dm(2)

    if(nproc_id_h/=ibox2)then
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

