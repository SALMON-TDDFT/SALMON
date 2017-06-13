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
!============================ Hartree potential (Solve Poisson equation)
SUBROUTINE Hartree_cg(trho,tVh)
use salmon_parallel, only: nproc_id_global, nproc_size_global, nproc_group_h
use salmon_communication, only: comm_is_root
use mpi, only: mpi_wtime, mpi_double_precision, mpi_sum
use scf_data
use new_world_sub
use sendrecvh_sub
use allocate_mat_sub
use deallocate_mat_sub

implicit none
real(8) :: trho(mg_sta(1):mg_end(1),    &
               mg_sta(2):mg_end(2),      &
               mg_sta(3):mg_end(3))
real(8) :: tVh(mg_sta(1):mg_end(1),    &
               mg_sta(2):mg_end(2),      &
               mg_sta(3):mg_end(3))

integer,parameter :: maxiter=1000
integer :: ix,iy,iz,iter
real(8) :: sum1,sum2,ak,ck
real(8) :: tottmp
real(8) :: totbox
real(8) :: wk2(ng_sta(1)-Ndh:ng_end(1)+Ndh,    &
               ng_sta(2)-Ndh:ng_end(2)+Ndh,      &
               ng_sta(3)-Ndh:ng_end(3)+Ndh)
real(8) :: rlap_wk(ng_sta(1):ng_end(1),    &
               ng_sta(2):ng_end(2),      &
               ng_sta(3):ng_end(3))
real(8) :: zk(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) 
real(8) :: tk(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) 
real(8) :: pk(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)) 
integer :: ierr

iwk_size=12
call make_iwksta_iwkend

call hartree_boundary(trho,wk2)

!------------------------- C-G minimization

!$OMP parallel do
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  wk2(ix,iy,iz)=tVh(ix,iy,iz)
  zk(ix,iy,iz) = -4.d0*Pi*trho(ix,iy,iz)
end do
end do
end do

call sendrecvh(wk2)

call calc_laplacianh(wk2,rlap_wk)

!$OMP parallel do
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  zk(ix,iy,iz)=zk(ix,iy,iz)-rlap_wk(ix,iy,iz)
end do
end do
end do

!$OMP parallel do
do iz=ng_sta(3)-Ndh,ng_end(3)+Ndh
do iy=ng_sta(2)-Ndh,ng_end(2)+Ndh
do ix=ng_sta(1)-Ndh,ng_end(1)+Ndh
  wk2(ix,iy,iz)=0.d0 
end do
end do
end do
!$OMP parallel do
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  pk(ix,iy,iz)=zk(ix,iy,iz)
end do
end do
end do

sum1=0.d0
!$OMP parallel do reduction(+ : sum1)
do iz=ng_sta(3),ng_end(3)
do iy=ng_sta(2),ng_end(2)
do ix=ng_sta(1),ng_end(1)
  sum1=sum1+zk(ix,iy,iz)**2*Hvol
end do
end do
end do

if(nproc_size_global==1)then
else
  elp3(201)=MPI_Wtime()
  call MPI_ALLREDUCE(sum1,sum2,1, &
                  MPI_DOUBLE_PRECISION,      &
                  MPI_SUM,nproc_group_h,IERR)
  elp3(202)=MPI_Wtime()
  elp3(254)=elp3(254)+elp3(202)-elp3(201)
  sum1=sum2
end if

Iteration : do iter=1,maxiter


!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    wk2(ix,iy,iz)=pk(ix,iy,iz)
  end do
  end do
  end do

  call sendrecvh(wk2)

  call calc_laplacianh(wk2,rlap_wk)

!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    tk(ix,iy,iz)=rlap_wk(ix,iy,iz)
  end do
  end do
  end do

  totbox=0d0
!$OMP parallel do reduction(+ : totbox)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    totbox=totbox+(zk(ix,iy,iz)*tk(ix,iy,iz))
  end do
  end do
  end do

  if(nproc_size_global==1)then
    tottmp=totbox
  else
    elp3(201)=MPI_Wtime()
    call MPI_Allreduce(totbox,tottmp,1,MPI_DOUBLE_PRECISION,      &
                  MPI_SUM,nproc_group_h,ierr)
    elp3(202)=MPI_Wtime()
    elp3(255)=elp3(255)+elp3(202)-elp3(201)
  end if

  ak=sum1/tottmp/Hvol

!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
     tVh(ix,iy,iz)=tVh(ix,iy,iz)+ak*pk(ix,iy,iz)
     zk(ix,iy,iz)=zk(ix,iy,iz)-ak*tk(ix,iy,iz)
  end do
  end do
  end do

  totbox=0d0
!$OMP parallel do reduction(+ : totbox)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    totbox=totbox+zk(ix,iy,iz)**2
  end do
  end do
  end do

  if(nproc_size_global==1)then
    tottmp=totbox
  else
    elp3(201)=MPI_Wtime()
    call MPI_Allreduce(totbox,tottmp,1,MPI_DOUBLE_PRECISION,      &
                  MPI_SUM,nproc_group_h,ierr)
    elp3(202)=MPI_Wtime()
    elp3(256)=elp3(256)+elp3(202)-elp3(201)
  end if

  sum2=tottmp*Hvol

  if ( abs(sum2) < Hconv*dble(lg_num(1)*lg_num(2)*lg_num(3)) ) exit

  ck=sum2/sum1 ; sum1=sum2

!$OMP parallel do
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    pk(ix,iy,iz)=zk(ix,iy,iz)+ck*pk(ix,iy,iz)
  end do
  end do
  end do
   
end do Iteration

iterVh=iter
if ( iterVh>maxiter .and. comm_is_root(nproc_id_global)) then
   write(*,*) "Warning:Vh iteration is not converged"
   write(*,'("||tVh(i)-tVh(i-1)||**2/(# of grids) = ",e15.8)') &
                              sum2/dble(lg_num(1)*lg_num(2)*lg_num(3))
end if

return

END SUBROUTINE Hartree_cg
