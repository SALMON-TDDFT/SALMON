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
MODULE nabla_sub

use scf_data
!use sendrecv_sub
implicit none 
INTERFACE calc_nabla

  MODULE PROCEDURE R_calc_nabla,C_calc_nabla

END INTERFACE

CONTAINS

!=======================================================================
SUBROUTINE R_calc_nabla(wk2,rnab_wkx,rnab_wky,rnab_wkz)
!$ use omp_lib

implicit none
integer :: ix,iy,iz
real(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))
real(8) :: rnab_wkx(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
real(8) :: rnab_wky(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
real(8) :: rnab_wkz(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))

!$OMP parallel do private(iz,iy,ix)
do iz=iwk3sta(3),iwk3end(3)
do iy=iwk3sta(2),iwk3end(2)
do ix=iwk3sta(1),iwk3end(1)
  rnab_wkx(ix,iy,iz)= 0.5d0*( wk2(ix+1,iy,iz)-wk2(ix-1,iy,iz) )
  rnab_wky(ix,iy,iz)= 0.5d0*( wk2(ix,iy+1,iz)-wk2(ix,iy-1,iz) )
  rnab_wkz(ix,iy,iz)= 0.5d0*( wk2(ix,iy,iz+1)-wk2(ix,iy,iz-1) )
end do
end do
end do

return

END SUBROUTINE R_calc_nabla

!=======================================================================
SUBROUTINE C_calc_nabla(wk2,cnab_wkx,cnab_wky,cnab_wkz)
!$ use omp_lib

implicit none
integer :: ix,iy,iz
complex(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))
complex(8) :: cnab_wkx(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
complex(8) :: cnab_wky(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
complex(8) :: cnab_wkz(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))

!$OMP parallel do private(iz,iy,ix) 
do iz=iwk3sta(3),iwk3end(3)
do iy=iwk3sta(2),iwk3end(2)
do ix=iwk3sta(1),iwk3end(1)
  cnab_wkx(ix,iy,iz)= 0.5d0*( wk2(ix+1,iy,iz)-wk2(ix-1,iy,iz) )
  cnab_wky(ix,iy,iz)= 0.5d0*( wk2(ix,iy+1,iz)-wk2(ix,iy-1,iz) )
  cnab_wkz(ix,iy,iz)= 0.5d0*( wk2(ix,iy,iz+1)-wk2(ix,iy,iz-1) )
end do
end do
end do

return

END SUBROUTINE C_calc_nabla

END MODULE nabla_sub
