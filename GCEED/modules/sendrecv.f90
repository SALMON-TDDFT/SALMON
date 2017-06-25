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
module sendrecv_sub

  use scf_data
  use new_world_sub
  use init_sendrecv_sub
  
  interface sendrecv
  
    module procedure R_sendrecv, C_sendrecv
  
  end interface 

contains

!==================================================================================================

subroutine R_sendrecv(tpsi)
  use salmon_parallel, only: nproc_group_orbital
  use salmon_communication, only: comm_proc_null, comm_isend, comm_irecv, comm_wait_all
  implicit none
  real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd, &
                  mg_sta(3)-Nd:mg_end(3)+Nd)
  integer :: ix,iy,iz
  integer :: iup,idw,jup,jdw,kup,kdw
  integer :: ireq(12)
  integer :: icomm
  
  iup=iup_array(1)
  idw=idw_array(1)
  jup=jup_array(1)
  jdw=jdw_array(1)
  kup=kup_array(1)
  kdw=kdw_array(1)
  
  icomm=nproc_group_orbital
  
  !send from idw to iup
  
  if(iup/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      srmatbox1_x_3d(ix,iy,iz)=tpsi(mg_end(1)-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1)
    end do
    end do
    end do
  end if
  ireq(1) = comm_isend(srmatbox1_x_3d,iup,3,icomm)
  ireq(2) = comm_irecv(srmatbox2_x_3d,idw,3,icomm)
  
  !send from iup to idw
  
  if(idw/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      srmatbox3_x_3d(ix,iy,iz)=tpsi(mg_sta(1)+ix-1,iy+mg_sta(2)-1,iz+mg_sta(3)-1)
    end do
    end do
    end do
  end if
  ireq(3) = comm_isend(srmatbox3_x_3d,idw,4,icomm)
  ireq(4) = comm_irecv(srmatbox4_x_3d,iup,4,icomm)
  
  !send from jdw to jup
  
  if(jup/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      srmatbox1_y_3d(ix,iy,iz)=tpsi(ix+mg_sta(1)-1,mg_end(2)-Nd+iy,iz+mg_sta(3)-1)
    end do
    end do
    end do
  end if
  ireq(5) = comm_isend(srmatbox1_y_3d,jup,5,icomm)
  ireq(6) = comm_irecv(srmatbox2_y_3d,jdw,5,icomm)
  
  !send from jup to jdw
  
  if(jdw/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      srmatbox3_y_3d(ix,iy,iz)=tpsi(ix+mg_sta(1)-1,mg_sta(2)+iy-1,iz+mg_sta(3)-1)
    end do
    end do
    end do
  end if
  ireq(7) = comm_isend(srmatbox3_y_3d,jdw,6,icomm)
  ireq(8) = comm_irecv(srmatbox4_y_3d,jup,6,icomm)
  
  !send from kdw to kup
  
  if(kup/=comm_proc_null)then
    do iz=1,Nd
  !$OMP parallel do private(iy,ix)
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      srmatbox1_z_3d(ix,iy,iz)=tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)-Nd+iz)
    end do
    end do
    end do
  end if
  ireq( 9) = comm_isend(srmatbox1_z_3d,kup,7,icomm)
  ireq(10) = comm_irecv(srmatbox2_z_3d,kdw,7,icomm)
  
  !send from kup to kdw
  
  if(kdw/=comm_proc_null)then
    do iz=1,Nd
  !$OMP parallel do private(iy,ix)
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      srmatbox3_z_3d(ix,iy,iz)=tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)+iz-1)
    end do
    end do
    end do
  end if
  ireq(11) = comm_isend(srmatbox3_z_3d,kdw,8,icomm)
  ireq(12) = comm_irecv(srmatbox4_z_3d,kup,8,icomm)
  
  
  call comm_wait_all(ireq(1:2))
  if(idw/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      tpsi(mg_sta(1)-1-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1)=srmatbox2_x_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(3:4))
  if(iup/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      tpsi(mg_end(1)+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1)=srmatbox4_x_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(5:6))
  if(jdw/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,mg_sta(2)-1-Nd+iy,iz+mg_sta(3)-1)=srmatbox2_y_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(7:8))
  if(jup/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,mg_end(2)+iy,iz+mg_sta(3)-1)=srmatbox4_y_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(9:10))
  if(kdw/=comm_proc_null)then
    do iz=1,Nd
  !$OMP parallel do private(iy,ix)
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)-1-Nd+iz)=srmatbox2_z_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(11:12))
  if(kup/=comm_proc_null)then
    do iz=1,Nd
  !$OMP parallel do private(iy,ix)
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)+iz)=srmatbox4_z_3d(ix,iy,iz)
    end do
    end do
    end do
  end if

end subroutine R_sendrecv

!==================================================================================================

subroutine C_sendrecv(tpsi)
  use salmon_parallel, only: nproc_group_orbital
  use salmon_communication, only: comm_proc_null, comm_isend, comm_irecv, comm_wait_all
  implicit none
  complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd, &
                     mg_sta(3)-Nd:mg_end(3)+Nd)
  integer :: ix,iy,iz
  integer :: iup,idw,jup,jdw,kup,kdw
  integer :: icomm
  integer :: ireq(12)
  
  iup=iup_array(1)
  idw=idw_array(1)
  jup=jup_array(1)
  jdw=jdw_array(1)
  kup=kup_array(1)
  kdw=kdw_array(1)
  
  icomm=nproc_group_orbital
  
  !send from idw to iup
  
  if(iup/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      scmatbox1_x_3d(ix,iy,iz)=tpsi(mg_end(1)-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1)
    end do
    end do
    end do
  end if
  ireq(1) = comm_isend(scmatbox1_x_3d,iup,3,icomm)
  ireq(2) = comm_irecv(scmatbox2_x_3d,idw,3,icomm)
  
  !send from iup to idw
  
  if(idw/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      scmatbox3_x_3d(ix,iy,iz)=tpsi(mg_sta(1)+ix-1,iy+mg_sta(2)-1,iz+mg_sta(3)-1)
    end do
    end do
    end do
  end if
  ireq(3) = comm_isend(scmatbox3_x_3d,idw,4,icomm)
  ireq(4) = comm_irecv(scmatbox4_x_3d,iup,4,icomm)
  
  !send from jdw to jup
  
  if(jup/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      scmatbox1_y_3d(ix,iy,iz)=tpsi(ix+mg_sta(1)-1,mg_end(2)-Nd+iy,iz+mg_sta(3)-1)
    end do
    end do
    end do
  end if
  ireq(5) = comm_isend(scmatbox1_y_3d,jup,5,icomm)
  ireq(6) = comm_irecv(scmatbox2_y_3d,jdw,5,icomm)
  
  !send from jup to jdw
  
  if(jdw/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      scmatbox3_y_3d(ix,iy,iz)=tpsi(ix+mg_sta(1)-1,mg_sta(2)+iy-1,iz+mg_sta(3)-1)
    end do
    end do
    end do
  end if
  ireq(7) = comm_isend(scmatbox3_y_3d,jdw,6,icomm)
  ireq(8) = comm_irecv(scmatbox4_y_3d,jup,6,icomm)
  
  !send from kdw to kup
  
  if(kup/=comm_proc_null)then
    do iz=1,Nd
  !$OMP parallel do private(iy,ix)
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      scmatbox1_z_3d(ix,iy,iz)=tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)-Nd+iz)
    end do
    end do
    end do
  end if
  ireq( 9) = comm_isend(scmatbox1_z_3d,kup,7,icomm)
  ireq(10) = comm_irecv(scmatbox2_z_3d,kdw,7,icomm)
  
  !send from kup to kdw
  
  if(kdw/=comm_proc_null)then
    do iz=1,Nd
  !$OMP parallel do private(iy,ix)
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      scmatbox3_z_3d(ix,iy,iz)=tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)+iz-1)
    end do
    end do
    end do
  end if
  ireq(11) = comm_isend(scmatbox3_z_3d,kdw,8,icomm)
  ireq(12) = comm_irecv(scmatbox4_z_3d,kup,8,icomm)
  
  
  call comm_wait_all(ireq(1:2))
  if(idw/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      tpsi(mg_sta(1)-1-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1)=scmatbox2_x_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(3:4))
  if(iup/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      tpsi(mg_end(1)+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1)=scmatbox4_x_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(5:6))
  if(jdw/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,mg_sta(2)-1-Nd+iy,iz+mg_sta(3)-1)=scmatbox2_y_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(7:8))
  if(jup/=comm_proc_null)then
  !$OMP parallel do private(iz,iy,ix)
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,mg_end(2)+iy,iz+mg_sta(3)-1)=scmatbox4_y_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(9:10))
  if(kdw/=comm_proc_null)then
    do iz=1,Nd
  !$OMP parallel do private(iy,ix)
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)-1-Nd+iz)=scmatbox2_z_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  call comm_wait_all(ireq(11:12))
  if(kup/=comm_proc_null)then
    do iz=1,Nd
  !$OMP parallel do private(iy,ix)
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)+iz)=scmatbox4_z_3d(ix,iy,iz)
    end do
    end do
    end do
  end if
  
end subroutine C_sendrecv

!==================================================================================================

end module sendrecv_sub
