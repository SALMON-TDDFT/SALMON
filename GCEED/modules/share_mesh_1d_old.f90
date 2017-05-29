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
MODULE share_mesh_1d_old_sub
use scf_data
INTERFACE share_mesh_1d_old
  MODULE PROCEDURE R_share_mesh_1d_old,C_share_mesh_1d_old
END INTERFACE

CONTAINS

subroutine R_share_mesh_1d_old(rbox_array,icoo,ista,iend)
implicit none

integer :: icoo
integer :: ista,iend
integer :: iup,idw,jup,jdw,kup,kdw
integer :: isd_addr,irv_addr,jsd_addr,jrv_addr,ksd_addr,krv_addr
integer :: ii
integer :: istatus(MPI_STATUS_SIZE)
real(8) :: rbox_array(ista:iend)

if(icoo==1)then

  iup=myrank+1
  if(mod(myrank,nproc_Mxin(1))==nproc_Mxin(1)-1) iup=myrank+1-nproc_Mxin(1)

  idw=myrank-1
  if(mod(myrank,nproc_Mxin(1))==0) idw=myrank-1+nproc_Mxin(1)

  do ii=0,nproc_Mxin(1)-2 
    isd_addr=myrank/nproc_Mxin(1)*nproc_Mxin(1)+mod(myrank-ii+nproc_Mxin(1),nproc_Mxin(1))
    irv_addr=myrank/nproc_Mxin(1)*nproc_Mxin(1)+mod(myrank-1-ii+nproc_Mxin(1),nproc_Mxin(1))

    call mpi_sendrecv(rbox_array(ista_Mxin_old(1,isd_addr)),inum_Mxin_old(1,isd_addr),   &
                    MPI_DOUBLE_PRECISION,iup,1,   &
                    rbox_array(ista_Mxin_old(1,irv_addr)),inum_Mxin_old(1,irv_addr),   &
                    MPI_DOUBLE_PRECISION,idw,1,   &
                    MPI_COMM_WORLD,istatus,ierr)
  end do


else if(icoo==2)then

  jup=myrank+nproc_Mxin(1)
  if(mod(myrank/nproc_Mxin(1),nproc_Mxin(2))==nproc_Mxin(2)-1) jup=myrank+(1-nproc_Mxin(2))*nproc_Mxin(1)

  jdw=myrank-nproc_Mxin(1)
  if(mod(myrank/nproc_Mxin(1),nproc_Mxin(2))==0) jdw=myrank+(nproc_Mxin(2)-1)*nproc_Mxin(1)
  
  do ii=0,nproc_Mxin(2)-2 
    jsd_addr=myrank/(nproc_Mxin(1)*nproc_Mxin(2))*(nproc_Mxin(1)*nproc_Mxin(2))  & 
             +mod(myrank/nproc_Mxin(1)-ii+nproc_Mxin(2),nproc_Mxin(2))*nproc_Mxin(1)
    jrv_addr=myrank/(nproc_Mxin(1)*nproc_Mxin(2))*(nproc_Mxin(1)*nproc_Mxin(2))  & 
             +mod(myrank/nproc_Mxin(1)-1-ii+nproc_Mxin(2),nproc_Mxin(2))*nproc_Mxin(1)

    call mpi_sendrecv(rbox_array(ista_Mxin_old(2,jsd_addr)),inum_Mxin_old(2,jsd_addr),   &
                    MPI_DOUBLE_PRECISION,jup,1,   &
                    rbox_array(ista_Mxin_old(2,jrv_addr)),inum_Mxin_old(2,jrv_addr),   &
                    MPI_DOUBLE_PRECISION,jdw,1,   &
                    MPI_COMM_WORLD,istatus,ierr)
  end do


else if(icoo==3)then

  kup=myrank+nproc_Mxin(1)*nproc_Mxin(2)
  if(mod(myrank/(nproc_Mxin(1)*nproc_Mxin(2)),nproc_Mxin(3))==nproc_Mxin(3)-1)  &
    kup=myrank+(1-nproc_Mxin(3))*nproc_Mxin(1)*nproc_Mxin(2)

  kdw=myrank-nproc_Mxin(1)*nproc_Mxin(2)
  if(mod(myrank/(nproc_Mxin(1)*nproc_Mxin(2)),nproc_Mxin(3))==0)   &
    kdw=myrank+(nproc_Mxin(3)-1)*nproc_Mxin(1)*nproc_Mxin(2)
  
  do ii=0,nproc_Mxin(3)-2 
    ksd_addr=myrank/nproc_Mxin_mul*nproc_Mxin_mul  & 
             +mod(myrank/(nproc_Mxin(1)*nproc_Mxin(2))-ii+nproc_Mxin(3),nproc_Mxin(3))   &
              *nproc_Mxin(1)*nproc_Mxin(2)
    krv_addr=myrank/nproc_Mxin_mul*nproc_Mxin_mul  & 
             +mod(myrank/(nproc_Mxin(1)*nproc_Mxin(2))-1-ii+nproc_Mxin(3),nproc_Mxin(3))   &
              *nproc_Mxin(1)*nproc_Mxin(2)

    call mpi_sendrecv(rbox_array(ista_Mxin_old(3,ksd_addr)),inum_Mxin_old(3,ksd_addr),   &
                    MPI_DOUBLE_PRECISION,kup,1,   &
                    rbox_array(ista_Mxin_old(3,krv_addr)),inum_Mxin_old(3,krv_addr),   &
                    MPI_DOUBLE_PRECISION,kdw,1,   &
                    MPI_COMM_WORLD,istatus,ierr)
  end do


end if

end subroutine R_share_mesh_1d_old

subroutine C_share_mesh_1d_old(cbox_array,icoo,ista,iend)
implicit none

integer :: icoo
integer :: ista,iend
integer :: iup,idw,jup,jdw,kup,kdw
integer :: isd_addr,irv_addr,jsd_addr,jrv_addr,ksd_addr,krv_addr
integer :: ii
integer :: istatus(MPI_STATUS_SIZE)
complex(8) :: cbox_array(ista:iend)

if(icoo==1)then

  iup=myrank+1
  if(mod(myrank,nproc_Mxin(1))==nproc_Mxin(1)-1) iup=myrank+1-nproc_Mxin(1)

  idw=myrank-1
  if(mod(myrank,nproc_Mxin(1))==0) idw=myrank-1+nproc_Mxin(1)

  do ii=0,nproc_Mxin(1)-2 
    isd_addr=myrank/nproc_Mxin(1)*nproc_Mxin(1)+mod(myrank-ii+nproc_Mxin(1),nproc_Mxin(1))
    irv_addr=myrank/nproc_Mxin(1)*nproc_Mxin(1)+mod(myrank-1-ii+nproc_Mxin(1),nproc_Mxin(1))

    call mpi_sendrecv(cbox_array(ista_Mxin_old(1,isd_addr)),inum_Mxin_old(1,isd_addr),   &
                    MPI_DOUBLE_COMPLEX,iup,1,   &
                    cbox_array(ista_Mxin_old(1,irv_addr)),inum_Mxin_old(1,irv_addr),   &
                    MPI_DOUBLE_COMPLEX,idw,1,   &
                    MPI_COMM_WORLD,istatus,ierr)
  end do


else if(icoo==2)then

  jup=myrank+nproc_Mxin(1)
  if(mod(myrank/nproc_Mxin(1),nproc_Mxin(2))==nproc_Mxin(2)-1) jup=myrank+(1-nproc_Mxin(2))*nproc_Mxin(1)

  jdw=myrank-nproc_Mxin(1)
  if(mod(myrank/nproc_Mxin(1),nproc_Mxin(2))==0) jdw=myrank+(nproc_Mxin(2)-1)*nproc_Mxin(1)
  
  do ii=0,nproc_Mxin(2)-2 
    jsd_addr=myrank/(nproc_Mxin(1)*nproc_Mxin(2))*(nproc_Mxin(1)*nproc_Mxin(2))  & 
             +mod(myrank/nproc_Mxin(1)-ii+nproc_Mxin(2),nproc_Mxin(2))*nproc_Mxin(1)
    jrv_addr=myrank/(nproc_Mxin(1)*nproc_Mxin(2))*(nproc_Mxin(1)*nproc_Mxin(2))  & 
             +mod(myrank/nproc_Mxin(1)-1-ii+nproc_Mxin(2),nproc_Mxin(2))*nproc_Mxin(1)

    call mpi_sendrecv(cbox_array(ista_Mxin_old(2,jsd_addr)),inum_Mxin_old(2,jsd_addr),   &
                    MPI_DOUBLE_COMPLEX,jup,1,   &
                    cbox_array(ista_Mxin_old(2,jrv_addr)),inum_Mxin_old(2,jrv_addr),   &
                    MPI_DOUBLE_COMPLEX,jdw,1,   &
                    MPI_COMM_WORLD,istatus,ierr)
  end do


else if(icoo==3)then

  kup=myrank+nproc_Mxin(1)*nproc_Mxin(2)
  if(mod(myrank/(nproc_Mxin(1)*nproc_Mxin(2)),nproc_Mxin(3))==nproc_Mxin(3)-1)  &
    kup=myrank+(1-nproc_Mxin(3))*nproc_Mxin(1)*nproc_Mxin(2)

  kdw=myrank-nproc_Mxin(1)*nproc_Mxin(2)
  if(mod(myrank/(nproc_Mxin(1)*nproc_Mxin(2)),nproc_Mxin(3))==0)   &
    kdw=myrank+(nproc_Mxin(3)-1)*nproc_Mxin(1)*nproc_Mxin(2)
  
  do ii=0,nproc_Mxin(3)-2 
    ksd_addr=myrank/nproc_Mxin_mul*nproc_Mxin_mul  & 
             +mod(myrank/(nproc_Mxin(1)*nproc_Mxin(2))-ii+nproc_Mxin(3),nproc_Mxin(3))   &
              *nproc_Mxin(1)*nproc_Mxin(2)
    krv_addr=myrank/nproc_Mxin_mul*nproc_Mxin_mul  & 
             +mod(myrank/(nproc_Mxin(1)*nproc_Mxin(2))-1-ii+nproc_Mxin(3),nproc_Mxin(3))   &
              *nproc_Mxin(1)*nproc_Mxin(2)

    call mpi_sendrecv(cbox_array(ista_Mxin_old(3,ksd_addr)),inum_Mxin_old(3,ksd_addr),   &
                    MPI_DOUBLE_COMPLEX,kup,1,   &
                    cbox_array(ista_Mxin_old(3,krv_addr)),inum_Mxin_old(3,krv_addr),   &
                    MPI_DOUBLE_COMPLEX,kdw,1,   &
                    MPI_COMM_WORLD,istatus,ierr)
  end do


end if

end subroutine C_share_mesh_1d_old

END MODULE share_mesh_1d_old_sub
