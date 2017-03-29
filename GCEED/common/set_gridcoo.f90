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

subroutine set_gridcoo
use scf_data
!$use omp_lib
implicit none
integer :: ix,iy,iz
integer :: itNd

itNd=max(Nd,Ndh)
allocate(gridcoo(minval(lg_sta(:))-itNd:maxval(lg_end(:))+itNd,3))

if(imesh_oddeven==1)then
!$OMP parallel do
  do ix=lg_sta(1)-itNd,lg_end(1)+itNd
    gridcoo(ix,1)=dble(ix)*Hgs(1)
  end do
!$OMP parallel do
  do iy=lg_sta(2)-itNd,lg_end(2)+itNd
    gridcoo(iy,2)=dble(iy)*Hgs(2)
  end do
!$OMP parallel do
  do iz=lg_sta(3)-itNd,lg_end(3)+itNd
    gridcoo(iz,3)=dble(iz)*Hgs(3)
  end do
else if(imesh_oddeven==2)then
!$OMP parallel do
  do ix=lg_sta(1)-itNd,lg_end(1)+itNd
    gridcoo(ix,1)=(dble(ix)-0.5d0)*Hgs(1)
  end do
!$OMP parallel do
  do iy=lg_sta(2)-itNd,lg_end(2)+itNd
    gridcoo(iy,2)=(dble(iy)-0.5d0)*Hgs(2)
  end do
!$OMP parallel do
  do iz=lg_sta(3)-itNd,lg_end(3)+itNd
    gridcoo(iz,3)=(dble(iz)-0.5d0)*Hgs(3)
  end do
end if

end subroutine set_gridcoo
