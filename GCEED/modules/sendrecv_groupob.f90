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

module sendrecv_groupob_sub

use scf_data
use new_world_sub
use init_sendrecv_sub
use allocate_sendrecv_groupob_sub

interface sendrecv_groupob

  module procedure R_sendrecv_groupob, C_sendrecv_groupob

end interface 

contains

!==================================================================================================

subroutine R_sendrecv_groupob(tpsi)
implicit none
real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
integer :: ix,iy,iz,iob
integer :: iup,idw,jup,jdw,kup,kdw
integer :: icomm
integer :: ireq1,istatus1(MPI_STATUS_SIZE)
integer :: ireq2,istatus2(MPI_STATUS_SIZE)
integer :: ireq3,istatus3(MPI_STATUS_SIZE)
integer :: ireq4,istatus4(MPI_STATUS_SIZE)
integer :: ireq5,istatus5(MPI_STATUS_SIZE)
integer :: ireq6,istatus6(MPI_STATUS_SIZE)
integer :: ireq7,istatus7(MPI_STATUS_SIZE)
integer :: ireq8,istatus8(MPI_STATUS_SIZE)
integer :: ireq9,istatus9(MPI_STATUS_SIZE)
integer :: ireq10,istatus10(MPI_STATUS_SIZE)
integer :: ireq11,istatus11(MPI_STATUS_SIZE)
integer :: ireq12,istatus12(MPI_STATUS_SIZE)
integer :: istatus(MPI_STATUS_SIZE)
logical :: flag

iup=iup_array(1)
idw=idw_array(1)
jup=jup_array(1)
jdw=jdw_array(1)
kup=kup_array(1)
kdw=kdw_array(1)

icomm=newworld_comm_orbital

!send from idw to iup

if(iup/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      srmatbox1_x(ix,iy,iz,iob,1)=tpsi(mg_end(1)-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(srmatbox1_x,Nd*mg_num(2)*mg_num(3)*iobnum,MPI_DOUBLE_PRECISION,iup,3,icomm,ireq1,ierr)
call mpi_irecv(srmatbox2_x,Nd*mg_num(2)*mg_num(3)*iobnum,MPI_DOUBLE_PRECISION,idw,3,icomm,ireq2,ierr)
call mpi_test(ireq1,flag,istatus1,ierr)
call mpi_test(ireq2,flag,istatus2,ierr)

!send from iup to idw

if(idw/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      srmatbox3_x(ix,iy,iz,iob,1)=tpsi(mg_sta(1)+ix-1,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(srmatbox3_x,Nd*mg_num(2)*mg_num(3)*iobnum,MPI_DOUBLE_PRECISION,idw,4,icomm,ireq3,ierr)
call mpi_irecv(srmatbox4_x,Nd*mg_num(2)*mg_num(3)*iobnum,MPI_DOUBLE_PRECISION,iup,4,icomm,ireq4,ierr)
call mpi_test(ireq3,flag,istatus3,ierr)
call mpi_test(ireq4,flag,istatus4,ierr)

!send from jdw to jup

if(jup/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      srmatbox1_y(ix,iy,iz,iob,1)=tpsi(ix+mg_sta(1)-1,mg_end(2)-Nd+iy,iz+mg_sta(3)-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(srmatbox1_y,mg_num(1)*Nd*mg_num(3)*iobnum,MPI_DOUBLE_PRECISION,jup,5,icomm,ireq5,ierr)
call mpi_irecv(srmatbox2_y,mg_num(1)*Nd*mg_num(3)*iobnum,MPI_DOUBLE_PRECISION,jdw,5,icomm,ireq6,ierr)
call mpi_test(ireq5,flag,istatus5,ierr)
call mpi_test(ireq6,flag,istatus6,ierr)

!send from jup to jdw

if(jdw/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      srmatbox3_y(ix,iy,iz,iob,1)=tpsi(ix+mg_sta(1)-1,mg_sta(2)+iy-1,iz+mg_sta(3)-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(srmatbox3_y,mg_num(1)*Nd*mg_num(3)*iobnum,MPI_DOUBLE_PRECISION,jdw,6,icomm,ireq7,ierr)
call mpi_irecv(srmatbox4_y,mg_num(1)*Nd*mg_num(3)*iobnum,MPI_DOUBLE_PRECISION,jup,6,icomm,ireq8,ierr)
call mpi_test(ireq7,flag,istatus7,ierr)
call mpi_test(ireq8,flag,istatus8,ierr)

!send from kdw to kup

if(kup/=MPI_PROC_NULL)then
  do iob=1,iobnum
    do iz=1,Nd
!$OMP parallel do
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      srmatbox1_z(ix,iy,iz,iob,1)=tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)-Nd+iz,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(srmatbox1_z,mg_num(1)*mg_num(2)*Nd*iobnum,MPI_DOUBLE_PRECISION,kup,7,icomm,ireq9,ierr)
call mpi_irecv(srmatbox2_z,mg_num(1)*mg_num(2)*Nd*iobnum,MPI_DOUBLE_PRECISION,kdw,7,icomm,ireq10,ierr)
call mpi_test(ireq9,flag,istatus9,ierr)
call mpi_test(ireq10,flag,istatus10,ierr)

!send from kup to kdw

if(kdw/=MPI_PROC_NULL)then
  do iob=1,iobnum
    do iz=1,Nd
!$OMP parallel do
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      srmatbox3_z(ix,iy,iz,iob,1)=tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)+iz-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(srmatbox3_z,mg_num(1)*mg_num(2)*Nd*iobnum,MPI_DOUBLE_PRECISION,kdw,8,icomm,ireq11,ierr)
call mpi_irecv(srmatbox4_z,mg_num(1)*mg_num(2)*Nd*iobnum,MPI_DOUBLE_PRECISION,kup,8,icomm,ireq12,ierr)
call mpi_test(ireq11,flag,istatus11,ierr)
call mpi_test(ireq12,flag,istatus12,ierr)


call mpi_wait(ireq1,istatus1,ierr)
call mpi_wait(ireq2,istatus2,ierr)
if(idw/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      tpsi(mg_sta(1)-1-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)=srmatbox2_x(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq3,istatus3,ierr)
call mpi_wait(ireq4,istatus4,ierr)
if(iup/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      tpsi(mg_end(1)+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)=srmatbox4_x(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq5,istatus5,ierr)
call mpi_wait(ireq6,istatus6,ierr)
if(jdw/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,mg_sta(2)-1-Nd+iy,iz+mg_sta(3)-1,iob,1)=srmatbox2_y(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq7,istatus7,ierr)
call mpi_wait(ireq8,istatus8,ierr)
if(jup/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,mg_end(2)+iy,iz+mg_sta(3)-1,iob,1)=srmatbox4_y(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq9,istatus9,ierr)
call mpi_wait(ireq10,istatus10,ierr)
if(kdw/=MPI_PROC_NULL)then
  do iob=1,iobnum
    do iz=1,Nd
!$OMP parallel do
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)-1-Nd+iz,iob,1)=srmatbox2_z(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq11,istatus11,ierr)
call mpi_wait(ireq12,istatus12,ierr)
if(kup/=MPI_PROC_NULL)then
  do iob=1,iobnum
    do iz=1,Nd
!$OMP parallel do
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)+iz,iob,1)=srmatbox4_z(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

end subroutine R_sendrecv_groupob

!==================================================================================================

subroutine C_sendrecv_groupob(tpsi)
implicit none
complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                   mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
integer :: ix,iy,iz,iob
integer :: iup,idw,jup,jdw,kup,kdw
integer :: icomm
integer :: ireq1,istatus1(MPI_STATUS_SIZE)
integer :: ireq2,istatus2(MPI_STATUS_SIZE)
integer :: ireq3,istatus3(MPI_STATUS_SIZE)
integer :: ireq4,istatus4(MPI_STATUS_SIZE)
integer :: ireq5,istatus5(MPI_STATUS_SIZE)
integer :: ireq6,istatus6(MPI_STATUS_SIZE)
integer :: ireq7,istatus7(MPI_STATUS_SIZE)
integer :: ireq8,istatus8(MPI_STATUS_SIZE)
integer :: ireq9,istatus9(MPI_STATUS_SIZE)
integer :: ireq10,istatus10(MPI_STATUS_SIZE)
integer :: ireq11,istatus11(MPI_STATUS_SIZE)
integer :: ireq12,istatus12(MPI_STATUS_SIZE)
integer :: istatus(MPI_STATUS_SIZE)
logical :: flag

iup=iup_array(1)
idw=idw_array(1)
jup=jup_array(1)
jdw=jdw_array(1)
kup=kup_array(1)
kdw=kdw_array(1)

icomm=newworld_comm_orbital

!send from idw to iup

if(iup/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      scmatbox1_x(ix,iy,iz,iob,1)=tpsi(mg_end(1)-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(scmatbox1_x,Nd*mg_num(2)*mg_num(3)*iobnum,MPI_DOUBLE_COMPLEX,iup,3,icomm,ireq1,ierr)
call mpi_irecv(scmatbox2_x,Nd*mg_num(2)*mg_num(3)*iobnum,MPI_DOUBLE_COMPLEX,idw,3,icomm,ireq2,ierr)
call mpi_test(ireq1,flag,istatus1,ierr)
call mpi_test(ireq2,flag,istatus2,ierr)

!send from iup to idw

if(idw/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      scmatbox3_x(ix,iy,iz,iob,1)=tpsi(mg_sta(1)+ix-1,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(scmatbox3_x,Nd*mg_num(2)*mg_num(3)*iobnum,MPI_DOUBLE_COMPLEX,idw,4,icomm,ireq3,ierr)
call mpi_irecv(scmatbox4_x,Nd*mg_num(2)*mg_num(3)*iobnum,MPI_DOUBLE_COMPLEX,iup,4,icomm,ireq4,ierr)
call mpi_test(ireq3,flag,istatus3,ierr)
call mpi_test(ireq4,flag,istatus4,ierr)

!send from jdw to jup

if(jup/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      scmatbox1_y(ix,iy,iz,iob,1)=tpsi(ix+mg_sta(1)-1,mg_end(2)-Nd+iy,iz+mg_sta(3)-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(scmatbox1_y,mg_num(1)*Nd*mg_num(3)*iobnum,MPI_DOUBLE_COMPLEX,jup,5,icomm,ireq5,ierr)
call mpi_irecv(scmatbox2_y,mg_num(1)*Nd*mg_num(3)*iobnum,MPI_DOUBLE_COMPLEX,jdw,5,icomm,ireq6,ierr)
call mpi_test(ireq5,flag,istatus5,ierr)
call mpi_test(ireq6,flag,istatus6,ierr)

!send from jup to jdw

if(jdw/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      scmatbox3_y(ix,iy,iz,iob,1)=tpsi(ix+mg_sta(1)-1,mg_sta(2)+iy-1,iz+mg_sta(3)-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(scmatbox3_y,mg_num(1)*Nd*mg_num(3)*iobnum,MPI_DOUBLE_COMPLEX,jdw,6,icomm,ireq7,ierr)
call mpi_irecv(scmatbox4_y,mg_num(1)*Nd*mg_num(3)*iobnum,MPI_DOUBLE_COMPLEX,jup,6,icomm,ireq8,ierr)
call mpi_test(ireq7,flag,istatus7,ierr)
call mpi_test(ireq8,flag,istatus8,ierr)

!send from kdw to kup

if(kup/=MPI_PROC_NULL)then
  do iob=1,iobnum
    do iz=1,Nd
!$OMP parallel do
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      scmatbox1_z(ix,iy,iz,iob,1)=tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)-Nd+iz,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(scmatbox1_z,mg_num(1)*mg_num(2)*Nd*iobnum,MPI_DOUBLE_COMPLEX,kup,7,icomm,ireq9,ierr)
call mpi_irecv(scmatbox2_z,mg_num(1)*mg_num(2)*Nd*iobnum,MPI_DOUBLE_COMPLEX,kdw,7,icomm,ireq10,ierr)
call mpi_test(ireq9,flag,istatus9,ierr)
call mpi_test(ireq10,flag,istatus10,ierr)

!send from kup to kdw

if(kdw/=MPI_PROC_NULL)then
  do iob=1,iobnum
    do iz=1,Nd
!$OMP parallel do
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      scmatbox3_z(ix,iy,iz,iob,1)=tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)+iz-1,iob,1)
    end do
    end do
    end do
  end do
end if
call mpi_isend(scmatbox3_z,mg_num(1)*mg_num(2)*Nd*iobnum,MPI_DOUBLE_COMPLEX,kdw,8,icomm,ireq11,ierr)
call mpi_irecv(scmatbox4_z,mg_num(1)*mg_num(2)*Nd*iobnum,MPI_DOUBLE_COMPLEX,kup,8,icomm,ireq12,ierr)
call mpi_test(ireq11,flag,istatus11,ierr)
call mpi_test(ireq12,flag,istatus12,ierr)


call mpi_wait(ireq1,istatus1,ierr)
call mpi_wait(ireq2,istatus2,ierr)
if(idw/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      tpsi(mg_sta(1)-1-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)=scmatbox2_x(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq3,istatus3,ierr)
call mpi_wait(ireq4,istatus4,ierr)
if(iup/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,mg_num(2)
    do ix=1,Nd
      tpsi(mg_end(1)+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)=scmatbox4_x(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq5,istatus5,ierr)
call mpi_wait(ireq6,istatus6,ierr)
if(jdw/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,mg_sta(2)-1-Nd+iy,iz+mg_sta(3)-1,iob,1)=scmatbox2_y(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq7,istatus7,ierr)
call mpi_wait(ireq8,istatus8,ierr)
if(jup/=MPI_PROC_NULL)then
  do iob=1,iobnum
!$OMP parallel do
    do iz=1,mg_num(3)
    do iy=1,Nd
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,mg_end(2)+iy,iz+mg_sta(3)-1,iob,1)=scmatbox4_y(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq9,istatus9,ierr)
call mpi_wait(ireq10,istatus10,ierr)
if(kdw/=MPI_PROC_NULL)then
  do iob=1,iobnum
    do iz=1,Nd
!$OMP parallel do
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)-1-Nd+iz,iob,1)=scmatbox2_z(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

call mpi_wait(ireq11,istatus11,ierr)
call mpi_wait(ireq12,istatus12,ierr)
if(kup/=MPI_PROC_NULL)then
  do iob=1,iobnum
    do iz=1,Nd
!$OMP parallel do
    do iy=1,mg_num(2)
    do ix=1,mg_num(1)
      tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)+iz,iob,1)=scmatbox4_z(ix,iy,iz,iob,1)
    end do
    end do
    end do
  end do
end if

end subroutine C_sendrecv_groupob

!==================================================================================================

end module sendrecv_groupob_sub
