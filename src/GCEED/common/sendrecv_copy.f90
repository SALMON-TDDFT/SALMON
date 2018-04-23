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
subroutine sendrecv_copy(tpsi)
!$ use omp_lib
use scf_data
implicit none
integer :: ix,iy,iz,iob,iik
real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,k_sta:k_end)

!$omp parallel default(none) &
!$omp          shared(k_sta,k_end,iobnum,mg_sta,mg_end,tpsi,psi) &
!$omp          private(iik,iob,iz,iy,ix)
do iik=k_sta,k_end
do iob=1,iobnum
!$omp do collapse(2)
  do iz=mg_sta(3)-Nd,mg_end(3)+Nd
  do iy=mg_sta(2)-Nd,mg_end(2)+Nd
  do ix=mg_sta(1)-Nd,mg_end(1)+Nd
    tpsi(ix,iy,iz,iob,iik)=0.d0
  end do
  end do
  end do
!$omp end do
!$omp do collapse(2)
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    tpsi(ix,iy,iz,iob,iik)=psi(ix,iy,iz,iob,iik)
  end do
  end do
  end do
!$omp end do
end do
end do
!$omp end parallel

end subroutine sendrecv_copy
