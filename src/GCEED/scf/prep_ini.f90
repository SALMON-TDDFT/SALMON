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
use salmon_communication, only: comm_is_root, comm_bcast
use scf_data
implicit none
integer :: imol,jj
real(8),parameter :: rtmp=1.d-10

if(comm_is_root(nproc_id_global))then
  open(60,file=file_ini)
  read(60,*) file_OUT_ini
  read(60,*) H_ini
  read(60,*) (rLsize_ini(jj),jj=1,3)
  read(60,*) num_mol
  read(60,*) rlatcon
end if

call comm_bcast(file_OUT_ini, nproc_group_global)
call comm_bcast(H_ini,        nproc_group_global)
call comm_bcast(rLsize_ini,   nproc_group_global)
call comm_bcast(num_mol,      nproc_group_global)
call comm_bcast(rlatcon,      nproc_group_global)

allocate(coo_mol_ini(3,num_mol))

if(comm_is_root(nproc_id_global))then
  do imol=1,num_mol
    read(60,*) (coo_mol_ini(jj,imol),jj=1,3)
  end do
  close(60)
end if

call comm_bcast(coo_mol_ini,nproc_group_global)

H_ini=H_ini/a_B
rlatcon=rlatcon/a_B
rLsize_ini(:)=rLsize_ini(:)/a_B

coo_mol_ini(:,:)=rlatcon*coo_mol_ini(:,:)

lg_end_ini(:)=int((rLsize_ini(:)+rtmp)/2.d0/Hgs(:))
do jj=1,3
  select case(imesh_oddeven(jj))
    case(1)
      lg_sta_ini(jj)=-(int((rLsize_ini(jj)+rtmp)/2.d0/Hgs(jj)))
    case(2)
      lg_sta_ini(jj)=-(int((rLsize_ini(jj)+rtmp)/2.d0/Hgs(jj)))+1
  end select
end do
lg_num_ini(:)=lg_end_ini(:)-lg_sta_ini(:)+1

end subroutine prep_ini
