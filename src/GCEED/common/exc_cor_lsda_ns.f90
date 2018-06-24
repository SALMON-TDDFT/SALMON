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
subroutine exc_cor_lsda_ns
  use salmon_parallel, only: nproc_group_h
  use salmon_communication, only: comm_summation
  use builtin_pz_sp, only: exc_cor_pz_sp
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  implicit none
  integer :: ix,iy,iz
  real(8) :: tot_exc

  call exc_cor_pz_sp(mg_num(1)*mg_num(2)*mg_num(3), rho_s, exc_m_tmp, eexc_m_tmp, vxc_s)

  tot_exc=0.d0 
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    tot_exc=tot_exc+eexc_m_tmp(ix,iy,iz)
  end do
  end do
  end do
  tot_exc=tot_exc*hvol
 
  call comm_summation(tot_exc,Exc,nproc_group_h)

  return
end subroutine exc_cor_lsda_ns
