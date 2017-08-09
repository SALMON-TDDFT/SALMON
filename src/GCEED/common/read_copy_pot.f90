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
subroutine read_copy_pot(potbox,matbox_read,ig_sta,ig_end)
use salmon_parallel, only: nproc_id_global, nproc_group_global
use salmon_communication, only: comm_is_root, comm_bcast
use scf_data
implicit none
integer :: ig_sta(3),ig_end(3)
real(8) :: potbox(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
real(8) :: matbox_read(ig_sta(1):ig_end(1),ig_sta(2):ig_end(2),ig_sta(3):ig_end(3))
integer :: ix,iy,iz

if(comm_is_root(nproc_id_global))then
  read(96) ((( matbox_read(ix,iy,iz),ix=ig_sta(1),ig_end(1))   &
                                    ,iy=ig_sta(2),ig_end(2))   &
                                    ,iz=ig_sta(3),ig_end(3))
end if

call comm_bcast(matbox_read,nproc_group_global)

!$OMP parallel do private(iz,iy,ix) 
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  potbox(ix,iy,iz)=matbox_read(ix,iy,iz)
end do
end do
end do

end subroutine read_copy_pot
