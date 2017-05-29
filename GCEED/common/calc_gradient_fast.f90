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
subroutine calc_gradient_fast(tpsi,grad_wk)
use scf_data
use sendrecv_groupob_sub
implicit none
real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
real(8) :: grad_wk(mg_sta(1):mg_end(1)+1,mg_sta(2):mg_end(2),mg_sta(3):mg_end(3), &
                   1:iobnum,1,3)
integer :: ix,iy,iz,iob

call sendrecv_copy(tpsi)
call sendrecv_groupob(tpsi)

if(Nd==4)then
  do iob=1,iobnum
!$OMP parallel private(iz)
    do iz=mg_sta(3),mg_end(3)
!$OMP do
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      grad_wk(ix,iy,iz,iob,1,1) =  &
        +bN1/Hgs(1)*( tpsi(ix+1,iy,iz,iob,1) - tpsi(ix-1,iy,iz,iob,1) )    &
        +bN2/Hgs(1)*( tpsi(ix+2,iy,iz,iob,1) - tpsi(ix-2,iy,iz,iob,1) )    &
        +bN3/Hgs(1)*( tpsi(ix+3,iy,iz,iob,1) - tpsi(ix-3,iy,iz,iob,1) )    &
        +bN4/Hgs(1)*( tpsi(ix+4,iy,iz,iob,1) - tpsi(ix-4,iy,iz,iob,1) )
      grad_wk(ix,iy,iz,iob,1,2) =  &
        +bN1/Hgs(2)*( tpsi(ix,iy+1,iz,iob,1) - tpsi(ix,iy-1,iz,iob,1) )    &
        +bN2/Hgs(2)*( tpsi(ix,iy+2,iz,iob,1) - tpsi(ix,iy-2,iz,iob,1) )    &
        +bN3/Hgs(2)*( tpsi(ix,iy+3,iz,iob,1) - tpsi(ix,iy-3,iz,iob,1) )    &
        +bN4/Hgs(2)*( tpsi(ix,iy+4,iz,iob,1) - tpsi(ix,iy-4,iz,iob,1) )
    end do
    end do
!$OMP end do nowait
!$OMP do
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      grad_wk(ix,iy,iz,iob,1,3) =  &
        +bN1/Hgs(3)*( tpsi(ix,iy,iz+1,iob,1) - tpsi(ix,iy,iz-1,iob,1) )    &
        +bN2/Hgs(3)*( tpsi(ix,iy,iz+2,iob,1) - tpsi(ix,iy,iz-2,iob,1) )    &
        +bN3/Hgs(3)*( tpsi(ix,iy,iz+3,iob,1) - tpsi(ix,iy,iz-3,iob,1) )    &
        +bN4/Hgs(3)*( tpsi(ix,iy,iz+4,iob,1) - tpsi(ix,iy,iz-4,iob,1) )
    end do
    end do
!$OMP end do nowait
    end do
!$OMP end parallel
  end do
end if

end subroutine calc_gradient_fast
