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
module sendrecv_self_sub

  use scf_data
  use new_world_sub
  use init_sendrecv_sub
  
  interface sendrecv_self
  
    module procedure R_sendrecv_self, C_sendrecv_self
  
  end interface 

contains

!==================================================================================================

subroutine R_sendrecv_self(tpsi,padding)
  implicit none
  integer :: padding
  real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+padding, &
                  mg_sta(2)-Nd:mg_end(2)+Nd, &
                  mg_sta(3)-Nd:mg_end(3)+Nd)
  integer :: ix,iy,iz

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,Nd
    tpsi(mg_sta(1)-1-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1) = &
    tpsi(mg_end(1)  -Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,Nd
    tpsi(mg_end(1)+ix,  iy+mg_sta(2)-1,iz+mg_sta(3)-1) = &
    tpsi(mg_sta(1)+ix-1,iy+mg_sta(2)-1,iz+mg_sta(3)-1)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,mg_num(3)
  do iy=1,Nd
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,mg_sta(2)-1-Nd+iy,iz+mg_sta(3)-1) = &
    tpsi(ix+mg_sta(1)-1,mg_end(2)  -Nd+iy,iz+mg_sta(3)-1)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,mg_num(3)
  do iy=1,Nd
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,mg_end(2)+iy,  iz+mg_sta(3)-1) = &
    tpsi(ix+mg_sta(1)-1,mg_sta(2)+iy-1,iz+mg_sta(3)-1)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,Nd
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)-1-Nd+iz) = &
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)  -Nd+iz)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,Nd
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)+iz  ) = &
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)+iz-1)
  end do
  end do
  end do
!$omp end parallel do
end subroutine R_sendrecv_self

!==================================================================================================

subroutine C_sendrecv_self(tpsi,padding)
  implicit none
  integer    :: padding
  complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+padding, &
                     mg_sta(2)-Nd:mg_end(2)+Nd, &
                     mg_sta(3)-Nd:mg_end(3)+Nd)
  integer :: ix,iy,iz

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,Nd
    tpsi(mg_sta(1)-1-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1) = &
    tpsi(mg_end(1)  -Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,Nd
    tpsi(mg_end(1)+ix,  iy+mg_sta(2)-1,iz+mg_sta(3)-1) = &
    tpsi(mg_sta(1)+ix-1,iy+mg_sta(2)-1,iz+mg_sta(3)-1)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,mg_num(3)
  do iy=1,Nd
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,mg_sta(2)-1-Nd+iy,iz+mg_sta(3)-1) = &
    tpsi(ix+mg_sta(1)-1,mg_end(2)  -Nd+iy,iz+mg_sta(3)-1)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,mg_num(3)
  do iy=1,Nd
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,mg_end(2)+iy,  iz+mg_sta(3)-1) = &
    tpsi(ix+mg_sta(1)-1,mg_sta(2)+iy-1,iz+mg_sta(3)-1)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,Nd
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)-1-Nd+iz) = &
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)  -Nd+iz)
  end do
  end do
  end do
!$omp end parallel do

!$omp parallel do default(none) private(iz,iy,ix) shared(mg_sta,mg_end,mg_num,tpsi)
  do iz=1,Nd
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)+iz  ) = &
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)+iz-1)
  end do
  end do
  end do
!$omp end parallel do
end subroutine C_sendrecv_self

!==================================================================================================

end module sendrecv_self_sub
