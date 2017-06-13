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
subroutine calcJxyz
use salmon_parallel, only: nproc_size_global, nproc_id_global, nproc_group_global
use salmon_communication, only: comm_is_root
use mpi, only: mpi_integer, mpi_sum
use scf_data
use allocate_psl_sub
implicit none
integer :: iatom,jj,j2,ix,iy,iz
integer :: ikoa
integer :: jshift
integer :: numj1(MI,0:nproc_size_global-1),numj2(MI,0:nproc_size_global-1)
real(8) :: rr
integer :: ierr

if(iSCFRT==1)then
  if(comm_is_root(nproc_id_global))then
    write(*,*) "max( Mps(iatom) ) = ",maxMps
  end if
end if

numj1=0
do iatom=1,MI
  ikoa=Kion(iatom) ; jj=0
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    rr=sqrt((gridcoo(ix,1)-Rion(1,iatom))**2      &
          +(gridcoo(iy,2)-Rion(2,iatom))**2      &
          +(gridcoo(iz,3)-Rion(3,iatom))**2)
    if ( rr < Rps(ikoa) ) then
      jj=jj+1 ; Jxyz_tmp1(1,jj,iatom)=ix ; Jxyz_tmp1(2,jj,iatom)=iy ; Jxyz_tmp1(3,jj,iatom)=iz 
    end if
  end do
  end do
  end do
  numj1(iatom,nproc_id_global)=jj
end do

call MPI_Allreduce(numj1,numj2,MI*nproc_size_global,MPI_INTEGER,MPI_SUM,nproc_group_global,ierr)

Jxyz_tmp2=0
do iatom=1,MI
  Mps(iatom)=sum(numj2(iatom,:))
  if(comm_is_root(nproc_id_global))then
    jshift=0
  else
    jshift=sum(numj2(iatom,0:nproc_id_global-1))
  end if
  do j2=1,numj2(iatom,nproc_id_global)
    Jxyz_tmp2(:,j2+jshift,iatom)=Jxyz_tmp1(:,j2,iatom)
  end do
end do

call MPI_Allreduce(Jxyz_tmp2,Jxyz,3*maxMps*MI,MPI_INTEGER,MPI_SUM,nproc_group_global,ierr)

Jxxyyzz=0

if(iSCFRT==1)then
  do iatom=1,MI
    if(comm_is_root(nproc_id_global))then
      write(*,*) "Mps =", Mps(iatom)
    end if
  end do
end if

if(iSCFRT==1)then
  if(comm_is_root(nproc_id_global))then
    write(*,*) "Mlmps =",Mlmps
  end if
end if

return

end subroutine calcJxyz
