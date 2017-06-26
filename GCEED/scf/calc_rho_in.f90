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
subroutine calc_rho_in
!$ use omp_lib
use scf_data
implicit none
integer :: ix,iy,iz,is

!$OMP parallel do private(iz,iy,ix)
do iz=mg_sta(3),mg_end(3)
do iy=mg_sta(2),mg_end(2)
do ix=mg_sta(1),mg_end(1)
  rho_in(ix,iy,iz,num_rho_stock+1)=rho(ix,iy,iz)
end do
end do
end do

if(ilsda==1)then
  do is=1,2
!$OMP parallel do private(iz,iy,ix)
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      rho_s_in(ix,iy,iz,is,num_rho_stock+1)=rho_s(ix,iy,iz,is)
    end do
    end do
    end do
  end do
end if

end subroutine calc_rho_in
