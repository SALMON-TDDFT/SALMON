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
MODULE copy_psi_mesh_sub

use scf_data
use new_world_sub
integer :: icopy_psi_mesh

INTERFACE copy_psi_mesh
  MODULE PROCEDURE R_copy_psi_mesh, C_copy_psi_mesh 
END INTERFACE 

CONTAINS

!======================================================================
subroutine R_copy_psi_mesh(tpsi_mesh)
use salmon_parallel, only: nproc_group_global, nproc_group_h
use salmon_communication, only: comm_summation
implicit none

real(8) :: tpsi_mesh(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),   &
                     1:itotMST,1)
integer :: ix,iy,iz,iob,iik
integer :: iob_myob
integer :: icheck_corrkob
real(8),allocatable :: matbox(:,:,:)
real(8),allocatable :: matbox2(:,:,:)

allocate(matbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
allocate(matbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))

if(icopy_psi_mesh==1)then
  do iik=1,num_kpoints_rd
  do iob=1,itotMST
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,iik,icheck_corrkob)
    matbox=0.d0
    if(icheck_corrkob==1)then
!$OMP parallel do private(iz,iy,ix) 
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        matbox(ix,iy,iz)=psi(ix,iy,iz,iob_myob,iik)
      end do
      end do
      end do
    end if

    call comm_summation(matbox,matbox2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)

!$OMP parallel do private(iz,iy,ix) 
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      tpsi_mesh(ix,iy,iz,iob,iik)=matbox2(ix,iy,iz)
    end do
    end do
    end do
  end do
  end do
else if(icopy_psi_mesh==2)then
  do iik=1,num_kpoints_rd
  do iob=1,itotMST
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,iik,icheck_corrkob)
    matbox=0.d0
!$OMP parallel do private(iz,iy,ix) 
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      matbox(ix,iy,iz)=psi_mesh(ix,iy,iz,iob,iik)
    end do
    end do
    end do

    call comm_summation(matbox,matbox2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

    if(icheck_corrkob==1)then
!$OMP parallel do private(iz,iy,ix) 
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        psi(ix,iy,iz,iob_myob,iik)=matbox2(ix,iy,iz)
      end do
      end do
      end do
    end if
  end do
  end do
end if

deallocate(matbox,matbox2)

return

end subroutine R_copy_psi_mesh

!======================================================================
subroutine C_copy_psi_mesh(tpsi_mesh)
use salmon_parallel, only: nproc_group_global, nproc_group_h
use salmon_communication, only: comm_summation
implicit none

complex(8) :: tpsi_mesh(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3),   &
                        1:itotMST,1)
integer :: ix,iy,iz,iob,iob_myob,iik
integer :: icheck_corrkob
complex(8),allocatable :: matbox(:,:,:)
complex(8),allocatable :: matbox2(:,:,:)

allocate(matbox(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
allocate(matbox2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))

if(icopy_psi_mesh==1)then
  do iik=1,num_kpoints_rd
  do iob=1,itotMST
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,iik,icheck_corrkob)
    matbox=0.d0
    if(icheck_corrkob==1)then
!$OMP parallel do private(iz,iy,ix) 
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        matbox(ix,iy,iz)=zpsi(ix,iy,iz,iob_myob,iik)
      end do
      end do
      end do
    end if

    call comm_summation(matbox,matbox2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)

!$OMP parallel do private(iz,iy,ix) 
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      tpsi_mesh(ix,iy,iz,iob,iik)=matbox2(ix,iy,iz)
    end do
    end do
    end do
  end do
  end do
else if(icopy_psi_mesh==2)then
  do iik=1,num_kpoints_rd
  do iob=1,itotMST
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,iik,icheck_corrkob)
    matbox=0.d0
!$OMP parallel do private(iz,iy,ix) 
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      matbox(ix,iy,iz)=psi_mesh(ix,iy,iz,iob,iik)
    end do
    end do
    end do

    call comm_summation(matbox,matbox2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_h)

    if(icheck_corrkob==1)then
!$OMP parallel do private(iz,iy,ix)
      do iz=mg_sta(3),mg_end(3)
      do iy=mg_sta(2),mg_end(2)
      do ix=mg_sta(1),mg_end(1)
        zpsi(ix,iy,iz,iob_myob,iik)=matbox2(ix,iy,iz)
      end do
      end do
      end do
    end if
  end do
  end do
end if

deallocate(matbox,matbox2)

return

end subroutine C_copy_psi_mesh

!======================================================================
END MODULE copy_psi_mesh_sub
