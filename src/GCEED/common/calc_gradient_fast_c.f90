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
subroutine calc_gradient_fast_c(tzpsi,cgrad_wk)
use scf_data
use sendrecv_groupob_sub
implicit none
complex(8) :: tzpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                    mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,k_sta:k_end)
complex(8) :: cgrad_wk(mg_sta(1):mg_end(1)+1,mg_sta(2):mg_end(2),mg_sta(3):mg_end(3), &
                       1:iobnum,k_sta:k_end,3)
integer :: ix,iy,iz,iob,iik

call sendrecv_groupob(tzpsi)

if(Nd==4)then
  do iik=k_sta,k_end
  do iob=1,iobnum
!$OMP parallel private(iz)
    do iz=mg_sta(3),mg_end(3)
!$OMP do private(iy,ix) 
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      cgrad_wk(ix,iy,iz,iob,iik,1) =  &
        +bN1/Hgs(1)*( tzpsi(ix+1,iy,iz,iob,iik) - tzpsi(ix-1,iy,iz,iob,iik) )    &
        +bN2/Hgs(1)*( tzpsi(ix+2,iy,iz,iob,iik) - tzpsi(ix-2,iy,iz,iob,iik) )    &
        +bN3/Hgs(1)*( tzpsi(ix+3,iy,iz,iob,iik) - tzpsi(ix-3,iy,iz,iob,iik) )    &
        +bN4/Hgs(1)*( tzpsi(ix+4,iy,iz,iob,iik) - tzpsi(ix-4,iy,iz,iob,iik) )
      cgrad_wk(ix,iy,iz,iob,iik,2) =  &
        +bN1/Hgs(2)*( tzpsi(ix,iy+1,iz,iob,iik) - tzpsi(ix,iy-1,iz,iob,iik) )    &
        +bN2/Hgs(2)*( tzpsi(ix,iy+2,iz,iob,iik) - tzpsi(ix,iy-2,iz,iob,iik) )    &
        +bN3/Hgs(2)*( tzpsi(ix,iy+3,iz,iob,iik) - tzpsi(ix,iy-3,iz,iob,iik) )    &
        +bN4/Hgs(2)*( tzpsi(ix,iy+4,iz,iob,iik) - tzpsi(ix,iy-4,iz,iob,iik) )
    end do
    end do
!$OMP end do nowait
!$OMP do private(iy,ix) 
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      cgrad_wk(ix,iy,iz,iob,iik,3) =  &
        +bN1/Hgs(3)*( tzpsi(ix,iy,iz+1,iob,iik) - tzpsi(ix,iy,iz-1,iob,iik) )    &
        +bN2/Hgs(3)*( tzpsi(ix,iy,iz+2,iob,iik) - tzpsi(ix,iy,iz-2,iob,iik) )    &
        +bN3/Hgs(3)*( tzpsi(ix,iy,iz+3,iob,iik) - tzpsi(ix,iy,iz-3,iob,iik) )    &
        +bN4/Hgs(3)*( tzpsi(ix,iy,iz+4,iob,iik) - tzpsi(ix,iy,iz-4,iob,iik) )
    end do
    end do
!$OMP end do nowait
    end do
!$OMP end parallel
  end do
  end do
end if

end subroutine calc_gradient_fast_c
