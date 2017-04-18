! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine conv_core_exc_cor
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  implicit none
  integer :: ix,iy,iz
  character(20) :: xc_tmp
  integer :: ispin
  real(8) :: cval
  real(8) :: tot_exc

  xc_tmp='pz'
  ispin=0
  cval=0.d0

  do iz=1,ng_num(3)
  do iy=1,ng_num(2)
  do ix=1,ng_num(1)
    rho_tmp(ix,iy,iz)=rho(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1)
  end do
  end do
  end do

  call core_exc_cor(xc_tmp, ispin, cval, ng_num(1), ng_num(2), ng_num(3), Hvol,  &
                    rho_tmp, exc_dummy2, exc_dummy, exc_dummy3, exc_dummy3, &
                    exc_dummy, exc_dummy, vxc_tmp, exc_dummy2, exc_dummy2, tot_exc)

  do iz=1,ng_num(3)
  do iy=1,ng_num(2)
  do ix=1,ng_num(1)
    vxc(ng_sta(1)+ix-1,ng_sta(2)+iy-1,ng_sta(3)+iz-1)=vxc_tmp(ix,iy,iz)
  end do
  end do
  end do
 
  call MPI_ALLREDUCE(tot_exc,Exc,1,MPI_DOUBLE_PRECISION,      &
                     MPI_SUM,newworld_comm_h,IERR)

  return
end subroutine conv_core_exc_cor
