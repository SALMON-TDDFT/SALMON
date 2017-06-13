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
subroutine subdip(rNe,ifunc)
use salmon_parallel, only: nproc_group_h, nproc_id_global
use salmon_communication, only: comm_is_root
use mpi, only: mpi_double_precision, mpi_sum, mpi_wtime
use scf_data
use new_world_sub
use allocate_mat_sub
implicit none
integer :: ifunc
integer :: i1,ix,iy,iz,i2
real(8) :: rNe
real(8) :: rbox_array(10), rbox_arrayq(3, 3)
real(8) :: rbox_array2(10), rbox_arrayq2(3, 3)
real(8) :: rbox1, rbox1q, rbox1q12, rbox1q23, rbox1q31
real(8) :: fact
real(8) :: absr2
integer :: ierr
   
  elp3(526)=MPI_Wtime()

!$OMP parallel do
   do i1=1,4
     rbox_array(i1)=0.d0
   end do
   do i1=1,3
     rbox1=0.d0
!$OMP parallel do reduction( + : rbox1 )
     do iz=ng_sta(3),ng_end(3)
     do iy=ng_sta(2),ng_end(2)
     do ix=ng_sta(1),ng_end(1)
       rbox1=rbox1+vecR(i1,ix,iy,iz)*rho(ix,iy,iz)
     end do
     end do
     end do
     rbox_array(i1)=rbox1
   end do
   
   rbox1=0.d0
!$OMP parallel do reduction( + : rbox1 )
   do iz=ng_sta(3),ng_end(3)
   do iy=ng_sta(2),ng_end(2)
   do ix=ng_sta(1),ng_end(1)
     rbox1=rbox1+rho(ix,iy,iz)
   end do
   end do
   end do
   rbox_array(4)=rbox1
   
!-----------QUADRUPOLE-start----------

if(quadrupole=='y')then
   rho_diff(:,:,:) = rho(:,:,:)-rho0(:,:,:)

!$OMP parallel do
   do i1=1,3
     do i2=1,3
       rbox_arrayq(i1, i2)=0.d0
     end do
   end do

! diag-compornents

   do i1=1,3
     rbox1q=0.d0
!$OMP parallel do reduction( + : rbox1q )
     do iz=ng_sta(3),ng_end(3)
     do iy=ng_sta(2),ng_end(2)
     do ix=ng_sta(1),ng_end(1)
       absr2=vecR(1,ix,iy,iz)**2+vecR(2,ix,iy,iz)**2+vecR(3,ix,iy,iz)**2
       rbox1q=rbox1q+(3.d0*vecR(i1,ix,iy,iz)*vecR(i1,ix,iy,iz)-absr2)*rho_diff(ix,iy,iz)
     end do
     end do
     end do
     rbox_arrayq(i1, i1)=rbox1q
   end do

! non-diag compornents

     rbox1q12=0.d0
     rbox1q23=0.d0
     rbox1q31=0.d0
!$OMP parallel do reduction( + : rbox1q12,rbox1q23,rbox1q31 )
     do iz=ng_sta(3),ng_end(3)
     do iy=ng_sta(2),ng_end(2)
     do ix=ng_sta(1),ng_end(1)
       rbox1q12=rbox1q12+3.d0*vecR(1,ix,iy,iz)*vecR(2,ix,iy,iz)*rho_diff(ix,iy,iz)
       rbox1q23=rbox1q23+3.d0*vecR(2,ix,iy,iz)*vecR(3,ix,iy,iz)*rho_diff(ix,iy,iz)
       rbox1q31=rbox1q31+3.d0*vecR(3,ix,iy,iz)*vecR(1,ix,iy,iz)*rho_diff(ix,iy,iz)
     end do
     end do
     end do
     rbox_arrayq(1,2)=rbox1q12 ; rbox_arrayq(2,1)=rbox1q12
     rbox_arrayq(2,3)=rbox1q23 ; rbox_arrayq(3,2)=rbox1q23
     rbox_arrayq(3,1)=rbox1q31 ; rbox_arrayq(1,3)=rbox1q31

end if
   
   !-----------QUADRUPOLE-end----------

   elp3(761)=MPI_Wtime()
   call MPI_allreduce(rbox_array,rbox_array2,4,MPI_DOUBLE_PRECISION,MPI_SUM,      &
           nproc_group_h,ierr)
   call MPI_allreduce(rbox_arrayq,rbox_arrayq2,9,MPI_DOUBLE_PRECISION,MPI_SUM,      &
           nproc_group_h,ierr)
   elp3(762)=MPI_Wtime()
   elp3(784)=elp3(784)+elp3(762)-elp3(761)

   rNe=rbox_array2(4)*Hvol               ! Number of electrons
   if(ifunc==1)then
     Dp(1:3,itt)=rbox_array2(1:3)*Hgs(1:3)*Hvol-vecDs(1:3)
     do i1=1,3
       Qp(1:3,i1,itt)=rbox_arrayq2(1:3,i1)*Hgs(1:3)*Hvol
     end do
     rIe(itt)=rNe
   else if(ifunc==2)then
     Dp(1:3,itt-1)=rbox_array2(1:3)*Hgs(1:3)*Hvol-vecDs(1:3)
     do i1=1,3
       Qp(1:3,i1,itt-1)=rbox_arrayq2(1:3,i1)*Hgs(1:3)*Hvol
     end do
     rIe(itt-1)=rNe
   end if

  if(comm_is_root(nproc_id_global))then
    if(ifunc==1)then
      if(circular=='y')then
        write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8,i5,f16.8)')       &
          itt,dble(itt)*dt*2.41888d-2, (Dp(i1,itt)*a_B,i1=1,3), rNe, Etot*2d0*Ry,iterVh,dble(cumnum)
      else
        write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8,i5)')       &
          itt,dble(itt)*dt*2.41888d-2, (Dp(i1,itt)*a_B,i1=1,3), rNe, Etot*2d0*Ry,iterVh
      end if
    else if(ifunc==2)then
      if(circular=='y')then
        write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8,i5,f16.8)')       &
          itt-1,dble(itt-1)*dt*2.41888d-2, (Dp(i1,itt-1)*a_B,i1=1,3), rNe, Etot*2d0*Ry,iterVh,dble(cumnum)
      else
        write(*,'(i8,f14.8, 3e16.8, f15.8,f18.8,i5)')       &
          itt-1,dble(itt-1)*dt*2.41888d-2, (Dp(i1,itt-1)*a_B,i1=1,3), rNe, Etot*2d0*Ry,iterVh
      end if
    end if
  end if

  fact=1.d0

  if(ilsda==0)then
    if(rNe.lt.ifMST(1)*2.d0*10.d0*fact.and.rNe.gt.ifMST(1)*2.d0/10.d0*fact)then
      continue
    else
      write(*,*) nproc_id_global,"t=",itt
      write(*,*) nproc_id_global,"rbox1=",rbox1
      write(*,*) nproc_id_global,"Ne=",rNe
      stop
    end if
  else if(ilsda==1)then
    if(rNe.lt.(ifMST(1)+ifMST(2))*10.d0*fact.and.rNe.gt.(ifMST(1)+ifMST(2))/10.d0*fact)then
      continue
    else
      write(*,*) nproc_id_global,"t=",itt
      write(*,*) nproc_id_global,"rbox1=",rbox1
      write(*,*) nproc_id_global,"Ne=",rNe
      stop
    end if
  end if

  elp3(527)=MPI_Wtime()
  elp3(539)=elp3(539)+elp3(527)-elp3(526)

end subroutine subdip
