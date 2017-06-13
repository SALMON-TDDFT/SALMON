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
!=======================================================================
subroutine inner_product3(matbox1,matbox2,rbox2)
use salmon_parallel, only: nproc_group_orbital
use mpi, only: mpi_double_precision, mpi_sum, mpi_wtime
use scf_data
use new_world_sub
!$ use omp_lib
implicit none
integer :: ix,iy,iz
real(8) :: matbox1(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
real(8) :: matbox2(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
real(8) :: rbox,rbox2
integer :: ierr

rbox=0.d0
!$omp parallel do reduction(+ : rbox)
do iz=iwk3sta(3),iwk3end(3)
do iy=iwk3sta(2),iwk3end(2)
do ix=iwk3sta(1),iwk3end(1)
  rbox=rbox+matbox1(ix,iy,iz)*matbox2(ix,iy,iz)
end do
end do
end do

elp3(186)=MPI_Wtime()
call MPI_Allreduce(rbox,rbox2,1,MPI_DOUBLE_PRECISION,MPI_SUM,nproc_group_orbital,ierr)
elp3(187)=MPI_Wtime()
elp3(190)=elp3(190)+elp3(187)-elp3(186)

end subroutine inner_product3
