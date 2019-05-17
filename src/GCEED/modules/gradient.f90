!
!  Copyright 2017-2019 SALMON developers
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
MODULE gradient_sub

use scf_data
use gradient2_sub
use sendrecv_sub
implicit none 
INTERFACE calc_gradient

  MODULE PROCEDURE R_calc_gradient,C_calc_gradient

END INTERFACE

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE R_calc_gradient(wk,grad_wk)
!$ use omp_lib

implicit none
real(8) :: wk(iwksta(1):iwkend(1),iwksta(2):iwkend(2),iwksta(3):iwkend(3))
real(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))
real(8) :: grad_wk(3,iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))

integer :: ix,iy,iz

if(iwk_size==1.or.iwk_size==11)then
  wk2=0.d0
!$OMP parallel do private(iz,iy,ix) 
  do iz=iwksta(3),iwkend(3)
  do iy=iwksta(2),iwkend(2)
  do ix=iwksta(1),iwkend(1)
    wk2(ix,iy,iz)=wk(ix,iy,iz)
  end do
  end do
  end do

  call sendrecv(wk2)
  call calc_gradient2(wk2,grad_wk)

else if(iwk_size==2.or.iwk_size==12)then
  
  call sendrecv(wk)
  call calc_gradient2(wk,grad_wk)

end if

return

END SUBROUTINE R_calc_gradient

!=======================================================================

SUBROUTINE C_calc_gradient(wk,grad_wk)
!$ use omp_lib

implicit none
complex(8) :: wk(iwksta(1):iwkend(1),iwksta(2):iwkend(2),iwksta(3):iwkend(3))
complex(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))
complex(8) :: grad_wk(3,iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))

integer :: ix,iy,iz

if(iwk_size==1.or.iwk_size==11)then
  wk2=0.d0
!$OMP parallel do private(iz,iy,ix) 
  do iz=iwksta(3),iwkend(3)
  do iy=iwksta(2),iwkend(2)
  do ix=iwksta(1),iwkend(1)
    wk2(ix,iy,iz)=wk(ix,iy,iz)
  end do
  end do
  end do

  call sendrecv(wk2)
  call calc_gradient2(wk2,grad_wk)

else if(iwk_size==2.or.iwk_size==12)then
  
  call sendrecv(wk)
  call calc_gradient2(wk,grad_wk)

end if

return

END SUBROUTINE C_calc_gradient

END MODULE gradient_sub


