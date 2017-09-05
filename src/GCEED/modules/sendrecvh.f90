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
MODULE sendrecvh_sub

use scf_data
use init_sendrecv_sub
implicit none

INTERFACE sendrecvh

  module procedure R_sendrecvh
  ! C_sendrecvh is never used in GCEED part.

END INTERFACE

CONTAINS

!======================================================================
SUBROUTINE R_sendrecvh(wk2)
use salmon_communication, only: comm_proc_null, comm_start_all, comm_wait_all
use scf_data
use init_sendrecv_sub
use new_world_sub
use pack_unpack
use persistent_comm, only: nreqs_rh

implicit none
real(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))

integer :: iup,idw,jup,jdw,kup,kdw
integer :: ireqs(12)

if(iwk_size>=1.and.iwk_size<=3)then

  iup=iup_array(1)
  idw=idw_array(1)
  jup=jup_array(1)
  jdw=jdw_array(1)
  kup=kup_array(1)
  kdw=kdw_array(1)

  ireqs(:) = nreqs_rh(:,1)

else if(iwk_size>=11.and.iwk_size<=13)then

  iup=iup_array(2)
  idw=idw_array(2)
  jup=jup_array(2)
  jdw=jdw_array(2)
  kup=kup_array(2)
  kdw=kdw_array(2)

  ireqs(:) = nreqs_rh(:,2)

else if(iwk_size>=31.and.iwk_size<=33)then

  iup=iup_array(4)
  idw=idw_array(4)
  jup=jup_array(4)
  jdw=jdw_array(4)
  kup=kup_array(4)
  kdw=kdw_array(4)

  ireqs(:) = nreqs_rh(:,4)

end if

!=====================================================================-

!send from idw to iup
if(iup/=comm_proc_null)then
    call copy_data(wk2(iwk3end(1)-Ndh+1:iwk3end(1),               &
                       iwk3sta(2)      :iwk3sta(2)+iwk3num(2)-1,  &
                       iwk3sta(3)      :iwk3sta(3)+iwk3num(3)-1), &
                   rmatbox1_x_h)
end if
call comm_start_all(ireqs(1:2))

!send from iup to idw
if(idw/=comm_proc_null)then
    call copy_data(wk2(iwk3sta(1):iwk3sta(1)+Ndh-1,         &
                       iwk3sta(2):iwk3sta(2)+iwk3num(2)-1,  &
                       iwk3sta(3):iwk3sta(3)+iwk3num(3)-1), &
                   rmatbox3_x_h)
end if
call comm_start_all(ireqs(3:4))

!=====================================================================-

!send from jdw to jup
if(jup/=comm_proc_null)then
    call copy_data(wk2(iwk3sta(1)      :iwk3sta(1)+iwk3num(1)-1,  &
                       iwk3end(2)-Ndh+1:iwk3end(2),               &
                       iwk3sta(3)      :iwk3sta(3)+iwk3num(3)-1), &
                   rmatbox1_y_h)
end if
call comm_start_all(ireqs(5:6))

!send from jup to jdw
if(jdw/=comm_proc_null)then
    call copy_data(wk2(iwk3sta(1):iwk3sta(1)+iwk3num(1)-1,  &
                       iwk3sta(2):iwk3sta(2)+Ndh-1,         &
                       iwk3sta(3):iwk3sta(3)+iwk3num(3)-1), &
                   rmatbox3_y_h)
end if
call comm_start_all(ireqs(7:8))

!=====================================================================-

!send from kdw to kup
if(kup/=comm_proc_null)then
    call copy_data(wk2(iwk3sta(1)      :iwk3sta(1)+iwk3num(1)-1, &
                       iwk3sta(2)      :iwk3sta(2)+iwk3num(2)-1, &
                       iwk3end(3)-Ndh+1:iwk3end(3)),             &
                   rmatbox1_z_h)
end if
call comm_start_all(ireqs(9:10))

!send from kup to kdw
if(kdw/=comm_proc_null)then
    call copy_data(wk2(iwk3sta(1):iwk3sta(1)+iwk3num(1)-1,  &
                       iwk3sta(2):iwk3sta(2)+iwk3num(2)-1,  &
                       iwk3sta(3):iwk3sta(3)+Ndh-1),        &
                   rmatbox3_z_h)
end if
call comm_start_all(ireqs(11:12))

!=====================================================================-

!recv from idw to iup
call comm_wait_all(ireqs(1:2))
if(idw/=comm_proc_null)then
    call copy_data(rmatbox2_x_h,                               &
                   wk2(iwk3sta(1)-Ndh:iwk3sta(1),              &
                       iwk3sta(2)    :iwk3sta(2)+iwk3num(2)-1, &
                       iwk3sta(3)    :iwk3sta(3)+iwk3num(3)-1))
end if

!recv from iup to idw
call comm_wait_all(ireqs(3:4))
if(iup/=comm_proc_null)then
    call copy_data(rmatbox4_x_h,                             &
                   wk2(iwk3end(1)+1:iwk3end(1)+Ndh,          &
                       iwk3sta(2)  :iwk3sta(2)+iwk3num(2)-1, &
                       iwk3sta(3)  :iwk3sta(3)+iwk3num(3)-1))
end if

!=====================================================================-

!recv from jdw to jup
call comm_wait_all(ireqs(5:6))
if(jdw/=comm_proc_null)then
    call copy_data(rmatbox2_y_h,                               &
                   wk2(iwk3sta(1)    :iwk3sta(1)+iwk3num(1)-1, &
                       iwk3sta(2)-Ndh:iwk3sta(2),              &
                       iwk3sta(3)    :iwk3sta(3)+iwk3num(3)-1))
end if

!recv from jup to jdw
call comm_wait_all(ireqs(7:8))
if(jup/=comm_proc_null)then
    call copy_data(rmatbox4_y_h,                              &
                   wk2(iwk3sta(1)   :iwk3sta(1)+iwk3num(1)-1, &
                        iwk3end(2)+1:iwk3end(2)+Ndh,          &
                        iwk3sta(3)  :iwk3sta(3)+iwk3num(3)-1))
end if

!=====================================================================-

!recv from kdw to kup
call comm_wait_all(ireqs(9:10))
if(kdw/=comm_proc_null)then
    call copy_data(rmatbox2_z_h,                               &
                   wk2(iwk3sta(1)    :iwk3sta(1)+iwk3num(1)-1, &
                       iwk3sta(2)    :iwk3sta(2)+iwk3num(2)-1, &
                       iwk3sta(3)-Ndh:iwk3sta(3)))
end if

!recv from kup to kdw
call comm_wait_all(ireqs(11:12))
if(kup/=comm_proc_null)then
    call copy_data(rmatbox4_z_h,                             &
                   wk2(iwk3sta(1)  :iwk3sta(1)+iwk3num(1)-1, &
                       iwk3sta(2)  :iwk3sta(2)+iwk3num(2)-1, &
                       iwk3end(3)+1:iwk3end(3)+Ndh))
end if

return

END SUBROUTINE R_sendrecvh

!======================================================================
SUBROUTINE C_sendrecvh(wk2)
use salmon_parallel, only: nproc_group_global, nproc_group_h
use salmon_communication, only: comm_proc_null, comm_exchange
use scf_data
use init_sendrecv_sub
use new_world_sub

implicit none
complex(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))

integer :: ibox
integer :: ix,iy,iz
integer :: iup,idw,jup,jdw,kup,kdw
integer :: icomm

if(iwk_size>=1.and.iwk_size<=3)then

  iup=iup_array(1)
  idw=idw_array(1)
  jup=jup_array(1)
  jdw=jdw_array(1)
  kup=kup_array(1)
  kdw=kdw_array(1)

  icomm=nproc_group_global

else if(iwk_size>=11.and.iwk_size<=13)then

  iup=iup_array(2)
  idw=idw_array(2)
  jup=jup_array(2)
  jdw=jdw_array(2)
  kup=kup_array(2)
  kdw=kdw_array(2)

  icomm=nproc_group_h

else if(iwk_size>=31.and.iwk_size<=33)then

  iup=iup_array(4)
  idw=idw_array(4)
  jup=jup_array(4)
  jdw=jdw_array(4)
  kup=kup_array(4)
  kdw=kdw_array(4)

  icomm=nproc_group_h

end if

!send from idw to iup

if(iwk_size>=1.and.iwk_size<=3)then
  ibox=1
else if((iwk_size>=11.and.iwk_size<=13).or.(iwk_size>=31.and.iwk_size<=33))then
  ibox=13
end if
if(iup/=comm_proc_null)then
!$OMP parallel do private(iz,iy,ix)
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Ndh
    cmatbox1_x_h(ix,iy,iz)=wk2(iwk3end(1)-Ndh+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call comm_exchange(rmatbox1_x_h,iup,rmatbox2_x_h,idw,1,icomm)
if(idw/=comm_proc_null)then
!$OMP parallel do private(iz,iy,ix)
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Ndh
    wk2(iwk3sta(1)-1-Ndh+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)=cmatbox2_x_h(ix,iy,iz)
  end do
  end do
  end do
end if

!send from iup to idw

if(iwk_size>=1.and.iwk_size<=3)then
  ibox=3
else if((iwk_size>=11.and.iwk_size<=13).or.(iwk_size>=31.and.iwk_size<=33))then
  ibox=15
end if
if(idw/=comm_proc_null)then
!$OMP parallel do private(iz,iy,ix)
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Ndh
    cmatbox1_x_h(ix,iy,iz)=wk2(iwk3sta(1)+ix-1,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call comm_exchange(rmatbox1_x_h,idw,rmatbox2_x_h,iup,1,icomm)
if(iup/=comm_proc_null)then
!$OMP parallel do private(iz,iy,ix)
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Ndh
    wk2(iwk3end(1)+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)=cmatbox2_x_h(ix,iy,iz)
  end do
  end do
  end do
end if


!send from jdw to jup

if(iwk_size>=1.and.iwk_size<=3)then
  ibox=5
else if((iwk_size>=11.and.iwk_size<=13).or.(iwk_size>=31.and.iwk_size<=33))then
  ibox=17
end if
if(jup/=comm_proc_null)then
!$OMP parallel do private(iz,iy,ix)
  do iz=1,iwk3num(3)
  do iy=1,Ndh
  do ix=1,iwk3num(1)
    cmatbox1_y_h(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iwk3end(2)-Ndh+iy,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call comm_exchange(rmatbox1_y_h,jup,rmatbox2_y_h,jdw,1,icomm)
if(jdw/=comm_proc_null)then
!$OMP parallel do private(iz,iy,ix)
  do iz=1,iwk3num(3)
  do iy=1,Ndh
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iwk3sta(2)-1-Ndh+iy,iz+iwk3sta(3)-1)=cmatbox2_y_h(ix,iy,iz)
  end do
  end do
  end do
end if

!send from jup to jdw

if(iwk_size>=1.and.iwk_size<=3)then
  ibox=7
else if((iwk_size>=11.and.iwk_size<=13).or.(iwk_size>=31.and.iwk_size<=33))then
  ibox=19
end if
if(jdw/=comm_proc_null)then
!$OMP parallel do private(iz,iy,ix)
  do iz=1,iwk3num(3)
  do iy=1,Ndh
  do ix=1,iwk3num(1)
    cmatbox1_y_h(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iwk3sta(2)+iy-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call comm_exchange(rmatbox1_y_h,jdw,rmatbox2_y_h,jup,1,icomm)
if(jup/=comm_proc_null)then
!$OMP parallel do private(iz,iy,ix)
  do iz=1,iwk3num(3)
  do iy=1,Ndh
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iwk3end(2)+iy,iz+iwk3sta(3)-1)=cmatbox2_y_h(ix,iy,iz)
  end do
  end do
  end do
end if

!send from kdw to kup

if(iwk_size>=1.and.iwk_size<=3)then
  ibox=9
else if((iwk_size>=11.and.iwk_size<=13).or.(iwk_size>=31.and.iwk_size<=33))then
  ibox=21
end if
if(kup/=comm_proc_null)then
  do iz=1,Ndh
!$OMP parallel do private(iy,ix)
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    cmatbox1_z_h(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3end(3)-Ndh+iz)
  end do
  end do
  end do
end if
call comm_exchange(rmatbox1_z_h,kup,rmatbox2_z_h,kdw,1,icomm)
if(kdw/=comm_proc_null)then
  do iz=1,Ndh
!$OMP parallel do private(iy,ix)
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3sta(3)-1-Ndh+iz)=cmatbox2_z_h(ix,iy,iz)
  end do
  end do
  end do
end if

!send from kup to kdw

if(iwk_size>=1.and.iwk_size<=3)then
  ibox=11
else if((iwk_size>=11.and.iwk_size<=13).or.(iwk_size>=31.and.iwk_size<=33))then
  ibox=23
end if
if(kdw/=comm_proc_null)then
  do iz=1,Ndh
!$OMP parallel do private(iy,ix)
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    cmatbox1_z_h(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3sta(3)+iz-1)
  end do
  end do
  end do
end if
call comm_exchange(rmatbox1_z_h,kdw,rmatbox2_z_h,kup,1,icomm)
if(kup/=comm_proc_null)then
  do iz=1,Ndh
!$OMP parallel do private(iy,ix)
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3end(3)+iz)=cmatbox2_z_h(ix,iy,iz)
  end do
  end do
  end do
end if

return

END SUBROUTINE C_sendrecvh

END MODULE sendrecvh_sub
