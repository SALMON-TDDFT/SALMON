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
subroutine prep_ini
use salmon_parallel, only: nproc_id_global, nproc_group_global
use salmon_communication, only: comm_is_root
use mpi, only: mpi_character, mpi_double_precision, mpi_integer
use scf_data
implicit none
integer :: imol,jj
real(8),parameter :: epsilon=1.d-10
integer :: ierr

if(comm_is_root(nproc_id_global))then
  open(60,file=file_ini)
  read(60,*) file_OUT_ini
  read(60,*) H_ini
  read(60,*) (rLsize_ini(jj),jj=1,3)
  read(60,*) num_mol
  read(60,*) rlatcon
end if

call MPI_Bcast(file_OUT_ini,100,MPI_CHARACTER,0,nproc_group_global,ierr)
call MPI_Bcast(H_ini,1,MPI_DOUBLE_PRECISION,0,nproc_group_global,ierr)
call MPI_Bcast(rLsize_ini,3,MPI_DOUBLE_PRECISION,0,nproc_group_global,ierr)
call MPI_Bcast(num_mol,1,MPI_INTEGER,0,nproc_group_global,ierr)
call MPI_Bcast(rlatcon,1,MPI_DOUBLE_PRECISION,0,nproc_group_global,ierr)

allocate(coo_mol_ini(3,num_mol))

if(comm_is_root(nproc_id_global))then
  do imol=1,num_mol
    read(60,*) (coo_mol_ini(jj,imol),jj=1,3)
  end do
  close(60)
end if

call MPI_Bcast(coo_mol_ini,3*num_mol,MPI_DOUBLE_PRECISION,0,nproc_group_global,ierr)

H_ini=H_ini/a_B
rlatcon=rlatcon/a_B
rLsize_ini(:)=rLsize_ini(:)/a_B

coo_mol_ini(:,:)=rlatcon*coo_mol_ini(:,:)

lg_end_ini(:)=int((rLsize_ini(:)+epsilon)/2.d0/Hgs(:))
if(imesh_oddeven==1)then
  lg_sta_ini(:)=-(int((rLsize_ini(:)+epsilon)/2.d0/Hgs(:)))
else if(imesh_oddeven==2)then 
  lg_sta_ini(:)=-(int((rLsize_ini(:)+epsilon)/2.d0/Hgs(:)))+1
end if
lg_num_ini(:)=lg_end_ini(:)-lg_sta_ini(:)+1

end subroutine prep_ini
