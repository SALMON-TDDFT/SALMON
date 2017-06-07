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
subroutine writedns
  use scf_data
  use allocate_mat_sub
  implicit none
  integer :: ix,iy,iz
  character(30) :: suffix
  character(30) :: phys_quantity
  character(10) :: filenum

  !$OMP parallel do collapse(2)
  do iz=lg_sta(3),lg_end(3)
  do iy=lg_sta(2),lg_end(2)
  do ix=lg_sta(1),lg_end(1)
    matbox_l(ix,iy,iz)=0.d0
  end do
  end do
  end do

  !$OMP parallel do collapse(2)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    matbox_l(ix,iy,iz)=rho(ix,iy,iz)
  end do
  end do
  end do
  
  call MPI_Allreduce(matbox_l,matbox_l2, &
  &             lg_num(1)*lg_num(2)*lg_num(3), &
  &             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  if(iSCFRT==1)then 
    suffix = "dns"
  else if(iSCFRT==2)then
    write(filenum, '(i8)') itt
    suffix = "dns"//adjustl(filenum)
  end if
  phys_quantity = "dns"
  if(format3d=='avs')then
    call writeavs(103,suffix,matbox_l2)
  else if(format3d=='cube')then
    call writecube(103,suffix,phys_quantity,matbox_l2)
  end if

  if(iSCFRT==2)then
    !$OMP parallel do collapse(2)
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do

    !$OMP parallel do collapse(2)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      matbox_l(ix,iy,iz)=rho(ix,iy,iz)-rho0(ix,iy,iz)
    end do
    end do
    end do
  
    call MPI_Allreduce(matbox_l,matbox_l2, &
  &             lg_num(1)*lg_num(2)*lg_num(3), &
  &             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    write(filenum, '(i8)') itt
    suffix = "dnsdiff"//adjustl(filenum)
    phys_quantity = "dnsdiff"
    if(format3d=='avs')then
      call writeavs(103,suffix,matbox_l2)
    else if(format3d=='cube')then
      call writecube(103,suffix,phys_quantity,matbox_l2)
    end if
  end if
 
end subroutine writedns
