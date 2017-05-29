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
SUBROUTINE calc_laplacianh(wk2,lap_wk)
!$ use omp_lib
use scf_data

implicit none
integer :: ix,iy,iz,ist
real(8) :: wk2(iwk2sta(1):iwk2end(1),iwk2sta(2):iwk2end(2),iwk2sta(3):iwk2end(3))
real(8) :: lap_wk(iwk3sta(1):iwk3end(1),iwk3sta(2):iwk3end(2),iwk3sta(3):iwk3end(3))
real(8) :: f0
real(8) :: Hinv2(3)

f0=(1.d0/Hgs(1)**2   &
   +1.d0/Hgs(2)**2   &
   +1.d0/Hgs(3)**2)
Hinv2(1:3)=1.d0/Hgs(1:3)**2

!$OMP parallel do collapse(2) private(ist)
do iz=iwk3sta(3),iwk3end(3)
do iy=iwk3sta(2),iwk3end(2)
do ix=iwk3sta(1),iwk3end(1)
  lap_wk(ix,iy,iz)=cNmat(0,Nd)*f0*wk2(ix,iy,iz)
  do ist=1,Ndh
    lap_wk(ix,iy,iz)=lap_wk(ix,iy,iz)     &
           +cNmat(ist,Ndh)*( (wk2(ix+ist,iy,iz) + wk2(ix-ist,iy,iz))*Hinv2(1)       &
                            +(wk2(ix,iy+ist,iz) + wk2(ix,iy-ist,iz))*Hinv2(2)       &
                            +(wk2(ix,iy,iz+ist) + wk2(ix,iy,iz-ist))*Hinv2(3)  )
  end do
end do
end do
end do

return

END SUBROUTINE calc_laplacianh
