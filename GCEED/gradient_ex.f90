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

!=======================================================================
SUBROUTINE calc_gradient_ex(wk2,grad_wk)
!$ use omp_lib
use scf_data

implicit none
integer :: ix,iy,iz
complex(8) :: wk2(iwk2sta(1):iwk2end(1)+1,iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))
complex(8) :: grad_wk(3,iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))

if(Nd<=3)then
!$OMP parallel do
  do iz=iwk3sta(3),iwk3end(3)
  do iy=iwk3sta(2),iwk3end(2)
  do ix=iwk3sta(1),iwk3end(1)
    grad_wk(1,ix,iy,iz)=      &
      (wk2(ix+1,iy,iz)-wk2(ix-1,iy,iz))/2.d0/Hgs(1) 
    grad_wk(2,ix,iy,iz)=      &
      (wk2(ix,iy+1,iz)-wk2(ix,iy-1,iz))/2.d0/Hgs(2) 
    grad_wk(3,ix,iy,iz)=      &
      (wk2(ix,iy,iz+1)-wk2(ix,iy,iz-1))/2.d0/Hgs(3) 
  end do
  end do
  end do
else if(Nd==4)then
!$OMP parallel do
  do iz=iwk3sta(3),iwk3end(3)
  do iy=iwk3sta(2),iwk3end(2)
  do ix=iwk3sta(1),iwk3end(1)
    grad_wk(1,ix,iy,iz)=      &
      (bN1*(wk2(ix+1,iy,iz)-wk2(ix-1,iy,iz))      &
      +bN2*(wk2(ix+2,iy,iz)-wk2(ix-2,iy,iz))      &
      +bN3*(wk2(ix+3,iy,iz)-wk2(ix-3,iy,iz))      &
      +bN4*(wk2(ix+4,iy,iz)-wk2(ix-4,iy,iz)))/Hgs(1) 
    grad_wk(2,ix,iy,iz)=      &
      (bN1*(wk2(ix,iy+1,iz)-wk2(ix,iy-1,iz))      &
      +bN2*(wk2(ix,iy+2,iz)-wk2(ix,iy-2,iz))      &
      +bN3*(wk2(ix,iy+3,iz)-wk2(ix,iy-3,iz))      &
      +bN4*(wk2(ix,iy+4,iz)-wk2(ix,iy-4,iz)))/Hgs(2) 
    grad_wk(3,ix,iy,iz)=      &
      (bN1*(wk2(ix,iy,iz+1)-wk2(ix,iy,iz-1))      &
      +bN2*(wk2(ix,iy,iz+2)-wk2(ix,iy,iz-2))      &
      +bN3*(wk2(ix,iy,iz+3)-wk2(ix,iy,iz-3))      &
      +bN4*(wk2(ix,iy,iz+4)-wk2(ix,iy,iz-4)))/Hgs(3) 
  end do
  end do
  end do
end if

return

END SUBROUTINE calc_gradient_ex
