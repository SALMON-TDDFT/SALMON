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

  iup=nproc_id_global+1
  if(mod(nproc_id_global,nproc_Mxin(1))==nproc_Mxin(1)-1) iup=nproc_id_global+1-nproc_Mxin(1)

  idw=nproc_id_global-1
  if(mod(nproc_id_global,nproc_Mxin(1))==0) idw=nproc_id_global-1+nproc_Mxin(1)

  do ii=0,nproc_Mxin(1)-2 
    isd_addr=nproc_id_global/nproc_Mxin(1)*nproc_Mxin(1)+mod(nproc_id_global-ii+nproc_Mxin(1),nproc_Mxin(1))
    irv_addr=nproc_id_global/nproc_Mxin(1)*nproc_Mxin(1)+mod(nproc_id_global-1-ii+nproc_Mxin(1),nproc_Mxin(1))

    call mpi_sendrecv(rbox_array(ista_Mxin_old(1,isd_addr)),inum_Mxin_old(1,isd_addr),   &
                    MPI_DOUBLE_PRECISION,iup,1,   &
                    rbox_array(ista_Mxin_old(1,irv_addr)),inum_Mxin_old(1,irv_addr),   &
                    MPI_DOUBLE_PRECISION,idw,1,   &
                    nproc_group_global,istatus,ierr)
  end do


else if(icoo==2)then

  jup=nproc_id_global+nproc_Mxin(1)
  if(mod(nproc_id_global/nproc_Mxin(1),nproc_Mxin(2))==nproc_Mxin(2)-1) jup=nproc_id_global+(1-nproc_Mxin(2))*nproc_Mxin(1)

  jdw=nproc_id_global-nproc_Mxin(1)
  if(mod(nproc_id_global/nproc_Mxin(1),nproc_Mxin(2))==0) jdw=nproc_id_global+(nproc_Mxin(2)-1)*nproc_Mxin(1)
  
  do ii=0,nproc_Mxin(2)-2 
    jsd_addr=nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2))*(nproc_Mxin(1)*nproc_Mxin(2))  & 
             +mod(nproc_id_global/nproc_Mxin(1)-ii+nproc_Mxin(2),nproc_Mxin(2))*nproc_Mxin(1)
    jrv_addr=nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2))*(nproc_Mxin(1)*nproc_Mxin(2))  & 
             +mod(nproc_id_global/nproc_Mxin(1)-1-ii+nproc_Mxin(2),nproc_Mxin(2))*nproc_Mxin(1)

    call mpi_sendrecv(rbox_array(ista_Mxin_old(2,jsd_addr)),inum_Mxin_old(2,jsd_addr),   &
                    MPI_DOUBLE_PRECISION,jup,1,   &
                    rbox_array(ista_Mxin_old(2,jrv_addr)),inum_Mxin_old(2,jrv_addr),   &
                    MPI_DOUBLE_PRECISION,jdw,1,   &
                    nproc_group_global,istatus,ierr)
  end do


else if(icoo==3)then

  kup=nproc_id_global+nproc_Mxin(1)*nproc_Mxin(2)
  if(mod(nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2)),nproc_Mxin(3))==nproc_Mxin(3)-1)  &
    kup=nproc_id_global+(1-nproc_Mxin(3))*nproc_Mxin(1)*nproc_Mxin(2)

  kdw=nproc_id_global-nproc_Mxin(1)*nproc_Mxin(2)
  if(mod(nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2)),nproc_Mxin(3))==0)   &
    kdw=nproc_id_global+(nproc_Mxin(3)-1)*nproc_Mxin(1)*nproc_Mxin(2)
  
  do ii=0,nproc_Mxin(3)-2 
    ksd_addr=nproc_id_global/nproc_Mxin_mul*nproc_Mxin_mul  & 
             +mod(nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2))-ii+nproc_Mxin(3),nproc_Mxin(3))   &
              *nproc_Mxin(1)*nproc_Mxin(2)
    krv_addr=nproc_id_global/nproc_Mxin_mul*nproc_Mxin_mul  & 
             +mod(nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2))-1-ii+nproc_Mxin(3),nproc_Mxin(3))   &
              *nproc_Mxin(1)*nproc_Mxin(2)

    call mpi_sendrecv(rbox_array(ista_Mxin_old(3,ksd_addr)),inum_Mxin_old(3,ksd_addr),   &
                    MPI_DOUBLE_PRECISION,kup,1,   &
                    rbox_array(ista_Mxin_old(3,krv_addr)),inum_Mxin_old(3,krv_addr),   &
                    MPI_DOUBLE_PRECISION,kdw,1,   &
                    nproc_group_global,istatus,ierr)
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

  iup=nproc_id_global+1
  if(mod(nproc_id_global,nproc_Mxin(1))==nproc_Mxin(1)-1) iup=nproc_id_global+1-nproc_Mxin(1)

  idw=nproc_id_global-1
  if(mod(nproc_id_global,nproc_Mxin(1))==0) idw=nproc_id_global-1+nproc_Mxin(1)

  do ii=0,nproc_Mxin(1)-2 
    isd_addr=nproc_id_global/nproc_Mxin(1)*nproc_Mxin(1)+mod(nproc_id_global-ii+nproc_Mxin(1),nproc_Mxin(1))
    irv_addr=nproc_id_global/nproc_Mxin(1)*nproc_Mxin(1)+mod(nproc_id_global-1-ii+nproc_Mxin(1),nproc_Mxin(1))

    call mpi_sendrecv(cbox_array(ista_Mxin_old(1,isd_addr)),inum_Mxin_old(1,isd_addr),   &
                    MPI_DOUBLE_COMPLEX,iup,1,   &
                    cbox_array(ista_Mxin_old(1,irv_addr)),inum_Mxin_old(1,irv_addr),   &
                    MPI_DOUBLE_COMPLEX,idw,1,   &
                    nproc_group_global,istatus,ierr)
  end do


else if(icoo==2)then

  jup=nproc_id_global+nproc_Mxin(1)
  if(mod(nproc_id_global/nproc_Mxin(1),nproc_Mxin(2))==nproc_Mxin(2)-1) jup=nproc_id_global+(1-nproc_Mxin(2))*nproc_Mxin(1)

  jdw=nproc_id_global-nproc_Mxin(1)
  if(mod(nproc_id_global/nproc_Mxin(1),nproc_Mxin(2))==0) jdw=nproc_id_global+(nproc_Mxin(2)-1)*nproc_Mxin(1)
  
  do ii=0,nproc_Mxin(2)-2 
    jsd_addr=nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2))*(nproc_Mxin(1)*nproc_Mxin(2))  & 
             +mod(nproc_id_global/nproc_Mxin(1)-ii+nproc_Mxin(2),nproc_Mxin(2))*nproc_Mxin(1)
    jrv_addr=nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2))*(nproc_Mxin(1)*nproc_Mxin(2))  & 
             +mod(nproc_id_global/nproc_Mxin(1)-1-ii+nproc_Mxin(2),nproc_Mxin(2))*nproc_Mxin(1)

    call mpi_sendrecv(cbox_array(ista_Mxin_old(2,jsd_addr)),inum_Mxin_old(2,jsd_addr),   &
                    MPI_DOUBLE_COMPLEX,jup,1,   &
                    cbox_array(ista_Mxin_old(2,jrv_addr)),inum_Mxin_old(2,jrv_addr),   &
                    MPI_DOUBLE_COMPLEX,jdw,1,   &
                    nproc_group_global,istatus,ierr)
  end do


else if(icoo==3)then

  kup=nproc_id_global+nproc_Mxin(1)*nproc_Mxin(2)
  if(mod(nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2)),nproc_Mxin(3))==nproc_Mxin(3)-1)  &
    kup=nproc_id_global+(1-nproc_Mxin(3))*nproc_Mxin(1)*nproc_Mxin(2)

  kdw=nproc_id_global-nproc_Mxin(1)*nproc_Mxin(2)
  if(mod(nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2)),nproc_Mxin(3))==0)   &
    kdw=nproc_id_global+(nproc_Mxin(3)-1)*nproc_Mxin(1)*nproc_Mxin(2)
  
  do ii=0,nproc_Mxin(3)-2 
    ksd_addr=nproc_id_global/nproc_Mxin_mul*nproc_Mxin_mul  & 
             +mod(nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2))-ii+nproc_Mxin(3),nproc_Mxin(3))   &
              *nproc_Mxin(1)*nproc_Mxin(2)
    krv_addr=nproc_id_global/nproc_Mxin_mul*nproc_Mxin_mul  & 
             +mod(nproc_id_global/(nproc_Mxin(1)*nproc_Mxin(2))-1-ii+nproc_Mxin(3),nproc_Mxin(3))   &
              *nproc_Mxin(1)*nproc_Mxin(2)

    call mpi_sendrecv(cbox_array(ista_Mxin_old(3,ksd_addr)),inum_Mxin_old(3,ksd_addr),   &
                    MPI_DOUBLE_COMPLEX,kup,1,   &
                    cbox_array(ista_Mxin_old(3,krv_addr)),inum_Mxin_old(3,krv_addr),   &
                    MPI_DOUBLE_COMPLEX,kdw,1,   &
                    nproc_group_global,istatus,ierr)
  end do


end if

end subroutine C_share_mesh_1d_old

END MODULE share_mesh_1d_old_sub
