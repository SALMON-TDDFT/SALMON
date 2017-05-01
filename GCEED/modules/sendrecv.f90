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



MODULE sendrecv_sub

use scf_data
use init_sendrecv_sub
implicit none 
integer :: istatuses(MPI_STATUS_SIZE,12),ireqs(12)

INTERFACE sendrecv

   MODULE PROCEDURE R_sendrecv,C_sendrecv

END INTERFACE

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE R_sendrecv(wk2)
!$ use omp_lib
use new_world_sub

implicit none
real(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))

integer :: ibox
integer :: ix,iy,iz
integer :: iup,idw,jup,jdw,kup,kdw
integer :: istatus(MPI_STATUS_SIZE)
integer :: icomm


if(iwk_size>=1.and.iwk_size<=3)then

  iup=iup_array(1)
  idw=idw_array(1)
  jup=jup_array(1)
  jdw=jdw_array(1)
  kup=kup_array(1)
  kdw=kdw_array(1)

  icomm=newworld_comm_orbital

else if(iwk_size>=11.and.iwk_size<=13)then

  iup=iup_array(2)
  idw=idw_array(2)
  jup=jup_array(2)
  jdw=jdw_array(2)
  kup=kup_array(2)
  kdw=kdw_array(2)

  icomm=newworld_comm_h

else if(iwk_size>=31.and.iwk_size<=33)then

  iup=iup_array(4)
  idw=idw_array(4)
  jup=jup_array(4)
  jdw=jdw_array(4)
  kup=kup_array(4)
  kdw=kdw_array(4)

  icomm=newworld_comm_h

end if


!send from idw to iup

if(iwk_size>=1.and.iwk_size<=3)then
  ibox=1
else if((iwk_size>=11.and.iwk_size<=13).or.(iwk_size>=31.and.iwk_size<=33))then
  ibox=13
end if
if(iup/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Nd
    rmatbox1_x_s(ix,iy,iz)=wk2(iwk3end(1)-Nd+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(rmatbox1_x_s,Nd*iwk3num(2)*iwk3num(3),MPI_DOUBLE_PRECISION,iup,1,   &
                  rmatbox2_x_s,Nd*iwk3num(2)*iwk3num(3),MPI_DOUBLE_PRECISION,idw,1,icomm,istatus,ierr)
if(idw/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Nd
    wk2(iwk3sta(1)-1-Nd+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)=rmatbox2_x_s(ix,iy,iz)
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
if(idw/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Nd
    rmatbox1_x_s(ix,iy,iz)=wk2(iwk3sta(1)+ix-1,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(rmatbox1_x_s,Nd*iwk3num(2)*iwk3num(3),MPI_DOUBLE_PRECISION,idw,1,   &
                  rmatbox2_x_s,Nd*iwk3num(2)*iwk3num(3),MPI_DOUBLE_PRECISION,iup,1,icomm,istatus,ierr)
if(iup/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Nd
    wk2(iwk3end(1)+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)=rmatbox2_x_s(ix,iy,iz)
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
if(jup/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,Nd
  do ix=1,iwk3num(1)
    rmatbox1_y_s(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iwk3end(2)-Nd+iy,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(rmatbox1_y_s,iwk3num(1)*Nd*iwk3num(3),MPI_DOUBLE_PRECISION,jup,1,   &
                  rmatbox2_y_s,iwk3num(1)*Nd*iwk3num(3),MPI_DOUBLE_PRECISION,jdw,1,icomm,istatus,ierr)
if(jdw/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,Nd
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iwk3sta(2)-1-Nd+iy,iz+iwk3sta(3)-1)=rmatbox2_y_s(ix,iy,iz)
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
if(jdw/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,Nd
  do ix=1,iwk3num(1)
    rmatbox1_y_s(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iwk3sta(2)+iy-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(rmatbox1_y_s,iwk3num(1)*Nd*iwk3num(3),MPI_DOUBLE_PRECISION,jdw,1,   &
                  rmatbox2_y_s,iwk3num(1)*Nd*iwk3num(3),MPI_DOUBLE_PRECISION,jup,1,icomm,istatus,ierr)
if(jup/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,Nd
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iwk3end(2)+iy,iz+iwk3sta(3)-1)=rmatbox2_y_s(ix,iy,iz)
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
if(kup/=MPI_PROC_NULL)then
  do iz=1,Nd
!$OMP parallel do
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    rmatbox1_z_s(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3end(3)-Nd+iz)
  end do
  end do
  end do
end if
call mpi_sendrecv(rmatbox1_z_s,iwk3num(1)*iwk3num(2)*Nd,MPI_DOUBLE_PRECISION,kup,7,   &
                  rmatbox2_z_s,iwk3num(1)*iwk3num(2)*Nd,MPI_DOUBLE_PRECISION,kdw,7,icomm,istatus,ierr)
if(kdw/=MPI_PROC_NULL)then
  do iz=1,Nd
!$OMP parallel do
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3sta(3)-1-Nd+iz)=rmatbox2_z_s(ix,iy,iz)
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
if(kdw/=MPI_PROC_NULL)then
  do iz=1,Nd
!$OMP parallel do
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    rmatbox1_z_s(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3sta(3)+iz-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(rmatbox1_z_s,iwk3num(1)*iwk3num(2)*Nd,MPI_DOUBLE_PRECISION,kdw,1,   &
                  rmatbox2_z_s,iwk3num(1)*iwk3num(2)*Nd,MPI_DOUBLE_PRECISION,kup,1,icomm,istatus,ierr)
if(kup/=MPI_PROC_NULL)then
  do iz=1,Nd
!$OMP parallel do
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3end(3)+iz)=rmatbox2_z_s(ix,iy,iz)
  end do
  end do
  end do
end if

return

END SUBROUTINE R_sendrecv

!=======================================================================
!=======================================================================

SUBROUTINE C_sendrecv(wk2)
!$ use omp_lib
use new_world_sub

implicit none
complex(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))

integer :: ibox
integer :: ix,iy,iz
integer :: iup,idw,jup,jdw,kup,kdw
integer :: istatus(MPI_STATUS_SIZE)
integer :: icomm


if(iwk_size>=1.and.iwk_size<=3)then

  iup=iup_array(1)
  idw=idw_array(1)
  jup=jup_array(1)
  jdw=jdw_array(1)
  kup=kup_array(1)
  kdw=kdw_array(1)

  icomm=newworld_comm_orbital

else if(iwk_size>=11.and.iwk_size<=13)then

  iup=iup_array(2)
  idw=idw_array(2)
  jup=jup_array(2)
  jdw=jdw_array(2)
  kup=kup_array(2)
  kdw=kdw_array(2)

  icomm=newworld_comm_h

else if(iwk_size>=31.and.iwk_size<=33)then

  iup=iup_array(4)
  idw=idw_array(4)
  jup=jup_array(4)
  jdw=jdw_array(4)
  kup=kup_array(4)
  kdw=kdw_array(4)

  icomm=newworld_comm_h

end if


!send from idw to iup

if(iwk_size>=1.and.iwk_size<=3)then
  ibox=1
else if((iwk_size>=11.and.iwk_size<=13).or.(iwk_size>=31.and.iwk_size<=33))then
  ibox=13
end if
if(iup/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Nd
    cmatbox1_x_s(ix,iy,iz)=wk2(iwk3end(1)-Nd+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(cmatbox1_x_s,Nd*iwk3num(2)*iwk3num(3),MPI_DOUBLE_COMPLEX,iup,3,   &
                  cmatbox2_x_s,Nd*iwk3num(2)*iwk3num(3),MPI_DOUBLE_COMPLEX,idw,3,icomm,istatus,ierr)
if(idw/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Nd
    wk2(iwk3sta(1)-1-Nd+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)=cmatbox2_x_s(ix,iy,iz)
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
if(idw/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Nd
    cmatbox1_x_s(ix,iy,iz)=wk2(iwk3sta(1)+ix-1,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(cmatbox1_x_s,Nd*iwk3num(2)*iwk3num(3),MPI_DOUBLE_COMPLEX,idw,4,   &
                  cmatbox2_x_s,Nd*iwk3num(2)*iwk3num(3),MPI_DOUBLE_COMPLEX,iup,4,icomm,istatus,ierr)
if(iup/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,iwk3num(2)
  do ix=1,Nd
    wk2(iwk3end(1)+ix,iy+iwk3sta(2)-1,iz+iwk3sta(3)-1)=cmatbox2_x_s(ix,iy,iz)
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
if(jup/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,Nd
  do ix=1,iwk3num(1)
    cmatbox1_y_s(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iwk3end(2)-Nd+iy,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(cmatbox1_y_s,iwk3num(1)*Nd*iwk3num(3),MPI_DOUBLE_COMPLEX,jup,5,   &
                  cmatbox2_y_s,iwk3num(1)*Nd*iwk3num(3),MPI_DOUBLE_COMPLEX,jdw,5,icomm,istatus,ierr)
if(jdw/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,Nd
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iwk3sta(2)-1-Nd+iy,iz+iwk3sta(3)-1)=cmatbox2_y_s(ix,iy,iz)
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
if(jdw/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,Nd
  do ix=1,iwk3num(1)
    cmatbox1_y_s(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iwk3sta(2)+iy-1,iz+iwk3sta(3)-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(cmatbox1_y_s,iwk3num(1)*Nd*iwk3num(3),MPI_DOUBLE_COMPLEX,jdw,6,   &
                  cmatbox2_y_s,iwk3num(1)*Nd*iwk3num(3),MPI_DOUBLE_COMPLEX,jup,6,icomm,istatus,ierr)
if(jup/=MPI_PROC_NULL)then
!$OMP parallel do
  do iz=1,iwk3num(3)
  do iy=1,Nd
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iwk3end(2)+iy,iz+iwk3sta(3)-1)=cmatbox2_y_s(ix,iy,iz)
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
if(kup/=MPI_PROC_NULL)then
  do iz=1,Nd
!$OMP parallel do
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    cmatbox1_z_s(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3end(3)-Nd+iz)
  end do
  end do
  end do
end if
call mpi_sendrecv(cmatbox1_z_s,iwk3num(1)*iwk3num(2)*Nd,MPI_DOUBLE_COMPLEX,kup,7,   &
                  cmatbox2_z_s,iwk3num(1)*iwk3num(2)*Nd,MPI_DOUBLE_COMPLEX,kdw,7,icomm,istatus,ierr)
if(kdw/=MPI_PROC_NULL)then
  do iz=1,Nd
!$OMP parallel do
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3sta(3)-1-Nd+iz)=cmatbox2_z_s(ix,iy,iz)
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
if(kdw/=MPI_PROC_NULL)then
  do iz=1,Nd
!$OMP parallel do
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    cmatbox1_z_s(ix,iy,iz)=wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3sta(3)+iz-1)
  end do
  end do
  end do
end if
call mpi_sendrecv(cmatbox1_z_s,iwk3num(1)*iwk3num(2)*Nd,MPI_DOUBLE_COMPLEX,kdw,8,   &
                  cmatbox2_z_s,iwk3num(1)*iwk3num(2)*Nd,MPI_DOUBLE_COMPLEX,kup,8,icomm,istatus,ierr)
if(kup/=MPI_PROC_NULL)then
  do iz=1,Nd
!$OMP parallel do
  do iy=1,iwk3num(2)
  do ix=1,iwk3num(1)
    wk2(ix+iwk3sta(1)-1,iy+iwk3sta(2)-1,iwk3end(3)+iz)=cmatbox2_z_s(ix,iy,iz)
  end do
  end do
  end do
end if

return

END SUBROUTINE C_sendrecv

END MODULE sendrecv_sub

