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
MODULE calc_density_sub

use scf_data
use new_world_sub
use allocate_mat_sub
implicit none
INTERFACE calc_density

  MODULE PROCEDURE R_calc_density,C_calc_density

END INTERFACE

CONTAINS

subroutine R_calc_density(tpsi,ifunc)
use salmon_parallel, only: nproc_group_grid
use mpi, only: mpi_double_precision, mpi_sum
implicit none
real(8) :: tpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,1)
integer :: i1_allob
integer :: ix,iy,iz,i1,iss,iob
integer :: iob_myob,icorr_p
integer :: iob_start(2),iob_end(2)
integer :: ifunc
integer :: ierr

if(ilsda==0)then
  matbox_m=0d0
  do i1=1,iobnum
    call calc_allob(i1,i1_allob)
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      matbox_m(ix,iy,iz)=matbox_m(ix,iy,iz)+abs(tpsi(ix,iy,iz,i1,1))**2*rocc(i1_allob,1)*wtk(1)
    end do
    end do
    end do
  end do
 
  if(ifunc==1)then 
    call MPI_Allreduce(matbox_m,rho(mg_sta(1),mg_sta(2),mg_sta(3)),      &
               mg_num(1)*mg_num(2)*mg_num(3), &
               MPI_DOUBLE_PRECISION,MPI_SUM,      &
               nproc_group_grid,ierr)
  else if(ifunc==2)then
    call MPI_Allreduce(matbox_m,rho_out(mg_sta(1),mg_sta(2),mg_sta(3),num_rho_stock),      &
               mg_num(1)*mg_num(2)*mg_num(3), &
               MPI_DOUBLE_PRECISION,MPI_SUM,      &
               nproc_group_grid,ierr)
  end if

else if(ilsda==1)then
  iob_start(1)=1
  iob_end(1)=MST(1)
  iob_start(2)=MST(1)+1
  iob_end(2)=itotMST
  do iss=1,2
    matbox_m=0d0
    do iob=iob_start(iss),iob_end(iss)
      call calc_myob(iob,iob_myob)
      call check_corrkob(iob,icorr_p)
      if(icorr_p==1)then
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          matbox_m(ix,iy,iz)=matbox_m(ix,iy,iz)+abs(tpsi(ix,iy,iz,iob_myob,1))**2*rocc(iob,1)*wtk(1)
        end do
        end do
        end do
      end if
    end do
    if(ifunc==1)then
      call MPI_Allreduce(matbox_m,rho_s(mg_sta(1),mg_sta(2),mg_sta(3),iss),      &
                 mg_num(1)*mg_num(2)*mg_num(3), &
                 MPI_DOUBLE_PRECISION,MPI_SUM,      &
                 nproc_group_grid,ierr)
    else if(ifunc==2)then
      call MPI_Allreduce(matbox_m,rho_s_out(mg_sta(1),mg_sta(2),mg_sta(3),iss,num_rho_stock),      &
                 mg_num(1)*mg_num(2)*mg_num(3), &
                 MPI_DOUBLE_PRECISION,MPI_SUM,      &
                 nproc_group_grid,ierr)
    end if
  end do
  
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rho(ix,iy,iz)=rho_s(ix,iy,iz,1)+rho_s(ix,iy,iz,2)
  end do
  end do
  end do
end if

return

END SUBROUTINE R_calc_density


subroutine C_calc_density(tpsi,ifunc)
use salmon_parallel, only: nproc_group_grid
use mpi, only: mpi_double_precision, mpi_sum
implicit none
complex(8) :: tpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum,1)
integer :: i1_allob
integer :: ix,iy,iz,i1,iss,iob
integer :: iob_myob,icorr_p
integer :: iob_start(2),iob_end(2)
integer :: ifunc
integer :: ierr

if(ilsda==0)then
  matbox_m=0d0
  do i1=1,iobnum
    call calc_allob(i1,i1_allob)
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      matbox_m(ix,iy,iz)=matbox_m(ix,iy,iz)+abs(tpsi(ix,iy,iz,i1,1))**2*rocc(i1_allob,1)*wtk(1)
    end do
    end do
    end do
  end do
 
  if(ifunc==1)then 
    call MPI_Allreduce(matbox_m,rho(mg_sta(1),mg_sta(2),mg_sta(3)),      &
               mg_num(1)*mg_num(2)*mg_num(3), &
               MPI_DOUBLE_PRECISION,MPI_SUM,      &
               nproc_group_grid,ierr)
  else if(ifunc==2)then
    call MPI_Allreduce(matbox_m,rho_out(mg_sta(1),mg_sta(2),mg_sta(3),num_rho_stock),      &
               mg_num(1)*mg_num(2)*mg_num(3), &
               MPI_DOUBLE_PRECISION,MPI_SUM,      &
               nproc_group_grid,ierr)
  end if

else if(ilsda==1)then
  iob_start(1)=1
  iob_end(1)=MST(1)
  iob_start(2)=MST(1)+1
  iob_end(2)=itotMST
  do iss=1,2
    matbox_m=0d0
    do iob=iob_start(iss),iob_end(iss)
      call calc_myob(iob,iob_myob)
      call check_corrkob(iob,icorr_p)
      if(icorr_p==1)then
!$OMP parallel do
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          matbox_m(ix,iy,iz)=matbox_m(ix,iy,iz)+abs(tpsi(ix,iy,iz,iob_myob,1))**2*rocc(iob,1)*wtk(1)
        end do
        end do
        end do
      end if
    end do
    if(ifunc==1)then
      call MPI_Allreduce(matbox_m,rho_s(mg_sta(1),mg_sta(2),mg_sta(3),iss),      &
                 mg_num(1)*mg_num(2)*mg_num(3), &
                 MPI_DOUBLE_PRECISION,MPI_SUM,      &
                 nproc_group_grid,ierr)
    else if(ifunc==2)then
      call MPI_Allreduce(matbox_m,rho_s_out(mg_sta(1),mg_sta(2),mg_sta(3),iss,num_rho_stock),      &
                 mg_num(1)*mg_num(2)*mg_num(3), &
                 MPI_DOUBLE_PRECISION,MPI_SUM,      &
                 nproc_group_grid,ierr)
    end if
  end do
  
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rho(ix,iy,iz)=rho_s(ix,iy,iz,1)+rho_s(ix,iy,iz,2)
  end do
  end do
  end do
end if

return

END SUBROUTINE C_calc_density

END MODULE calc_density_sub
