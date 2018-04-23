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
subroutine calcJxyz_all_periodic
use salmon_parallel, only: nproc_group_global, nproc_id_global, nproc_size_global
use salmon_communication, only: comm_is_root, comm_summation
use scf_data
use read_pslfile_sub
use allocate_psl_sub
implicit none

integer :: ix,iy,iz,ak
integer :: jj,iatom
real(8) :: x,y,z

complex(8),parameter :: zI=(0.d0,1.d0)
real(8) :: r

integer :: iix,iiy,iiz

integer :: max_icell(3,MKI)

integer :: jshift
integer :: numj1(MI,0:nproc_size_global-1),numj2(MI,0:nproc_size_global-1)
integer :: j2

! nonlocal potential
if(comm_is_root(nproc_id_global))then
  write(*,*) ''
  write(*,*) '============nonlocal grid data=============='
endif

if(iperiodic==1)then
  do ak=1,MKI
    if(Rps(ak)<Hgs(1)*lg_num(1))then
      max_icell(1,ak)=2
    else
      max_icell(1,ak)=Rps(ak)/(Hgs(1)*lg_num(1))+2
    end if
    max_icell(2:3,ak)=0
  end do
else if(iperiodic==3)then
  do ak=1,MKI
    max_icell(1:3,ak)=Rps(ak)/(Hgs(1:3)*lg_num(1:3))+2
  end do
end if

numj1=0
do iatom=1,MI
  ak=Kion(iatom)
  jj=0
  do iiz=-max_icell(3,ak),max_icell(3,ak)
  do iiy=-max_icell(2,ak),max_icell(2,ak)
  do iix=-max_icell(1,ak),max_icell(1,ak)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      x=gridcoo(ix,1)+dble(iix)*lg_num(1)*Hgs(1)-Rion(1,iatom)
      y=gridcoo(iy,2)+dble(iiy)*lg_num(2)*Hgs(2)-Rion(2,iatom)
      z=gridcoo(iz,3)+dble(iiz)*lg_num(3)*Hgs(3)-Rion(3,iatom)
      r=sqrt(x*x+y*y+z*z)
      if (r<Rps(ak)) then
        jj=jj+1
        Jxyz_tmp1(1,jj,iatom)=ix
        Jxyz_tmp1(2,jj,iatom)=iy
        Jxyz_tmp1(3,jj,iatom)=iz
        Jxxyyzz_tmp1(1,jj,iatom)=iix
        Jxxyyzz_tmp1(2,jj,iatom)=iiy
        Jxxyyzz_tmp1(3,jj,iatom)=iiz
      endif
    enddo
    enddo
    enddo
  enddo
  enddo
  enddo
  numj1(iatom,nproc_id_global)=jj
end do

call comm_summation(numj1,numj2,MI*nproc_size_global,nproc_group_global)

Jxyz_tmp2=0
Jxxyyzz_tmp2=0
do iatom=1,MI
  Mps_all(iatom)=sum(numj2(iatom,:))
  if(comm_is_root(nproc_id_global))then
    jshift=0
  else
    jshift=sum(numj2(iatom,0:nproc_id_global-1))
  end if
  do j2=1,numj2(iatom,nproc_id_global)
    Jxyz_tmp2(:,j2+jshift,iatom)=Jxyz_tmp1(:,j2,iatom)
    Jxxyyzz_tmp2(:,j2+jshift,iatom)=Jxxyyzz_tmp1(:,j2,iatom)
  end do
end do
call comm_summation(Jxyz_tmp2,Jxyz_all,3*maxMps*MI,       nproc_group_global)
call comm_summation(Jxxyyzz_tmp2,Jxxyyzz_all,3*maxMps*MI, nproc_group_global)


if(iSCFRT==1)then
  do iatom=1,MI
    if(comm_is_root(nproc_id_global))then
      write(*,*) "Mps =", Mps_all(iatom)
    end if
  end do
end if

return

end subroutine calcJxyz_all_periodic
