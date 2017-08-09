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
module sendrecv_groupob_ngp_sub

use scf_data

interface sendrecv_groupob_ngp

  module procedure R_sendrecv_groupob_ngp, C_sendrecv_groupob_ngp

end interface 

contains

!==================================================================================================

subroutine R_sendrecv_groupob_ngp(tpsi)
implicit none
real(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
integer :: ix,iy,iz,iob

do iob=1,iobnum
!$OMP parallel do private(iz,iy,ix)
  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,Nd
    tpsi(mg_sta(1)-1-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)=    &
                              tpsi(mg_end(1)-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)
  end do
  end do
  end do
end do

!send from iup to idw

do iob=1,iobnum
!$OMP parallel do private(iz,iy,ix)
  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,Nd
    tpsi(mg_end(1)+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)=   &
                              tpsi(mg_sta(1)+ix-1,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)
  end do
  end do
  end do
end do

!send from jdw to jup

do iob=1,iobnum
!$OMP parallel do private(iz,iy,ix)
  do iz=1,mg_num(3)
  do iy=1,Nd
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,mg_sta(2)-1-Nd+iy,iz+mg_sta(3)-1,iob,1)=   &
                              tpsi(ix+mg_sta(1)-1,mg_end(2)-Nd+iy,iz+mg_sta(3)-1,iob,1)
  end do
  end do
  end do
end do

!send from jup to jdw

do iob=1,iobnum
!$OMP parallel do private(iz,iy,ix)
  do iz=1,mg_num(3)
  do iy=1,Nd
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,mg_end(2)+iy,iz+mg_sta(3)-1,iob,1)=   &
                              tpsi(ix+mg_sta(1)-1,mg_sta(2)+iy-1,iz+mg_sta(3)-1,iob,1)
  end do
  end do
  end do
end do

!send from kdw to kup

do iob=1,iobnum
  do iz=1,Nd
!$OMP parallel do private(iy,ix)
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)-1-Nd+iz,iob,1)=  &
                              tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)-Nd+iz,iob,1)
  end do
  end do
  end do
end do

!send from kup to kdw

do iob=1,iobnum
  do iz=1,Nd
!$OMP parallel do private(iy,ix)
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)+iz,iob,1)= &
                              tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)+iz-1,iob,1)
  end do
  end do
  end do
end do

end subroutine R_sendrecv_groupob_ngp

!==================================================================================================

subroutine C_sendrecv_groupob_ngp(tpsi)
implicit none
complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,mg_sta(2)-Nd:mg_end(2)+Nd, &
                   mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
integer :: ix,iy,iz,iob

do iob=1,iobnum
!$OMP parallel do private(iz, iy,ix)
  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,Nd
    tpsi(mg_sta(1)-1-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)=    &
                              tpsi(mg_end(1)-Nd+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)
  end do
  end do
  end do
end do

!send from iup to idw

do iob=1,iobnum
!$OMP parallel do private(iz, iy,ix)
  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,Nd
    tpsi(mg_end(1)+ix,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)=   &
                              tpsi(mg_sta(1)+ix-1,iy+mg_sta(2)-1,iz+mg_sta(3)-1,iob,1)
  end do
  end do
  end do
end do

!send from jdw to jup

do iob=1,iobnum
!$OMP parallel do private(iz, iy,ix)
  do iz=1,mg_num(3)
  do iy=1,Nd
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,mg_sta(2)-1-Nd+iy,iz+mg_sta(3)-1,iob,1)=   &
                              tpsi(ix+mg_sta(1)-1,mg_end(2)-Nd+iy,iz+mg_sta(3)-1,iob,1)
  end do
  end do
  end do
end do

!send from jup to jdw

do iob=1,iobnum
!$OMP parallel do private(iz, iy,ix)
  do iz=1,mg_num(3)
  do iy=1,Nd
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,mg_end(2)+iy,iz+mg_sta(3)-1,iob,1)=   &
                              tpsi(ix+mg_sta(1)-1,mg_sta(2)+iy-1,iz+mg_sta(3)-1,iob,1)
  end do
  end do
  end do
end do

!send from kdw to kup

do iob=1,iobnum
  do iz=1,Nd
!$OMP parallel do private(iy,ix)
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)-1-Nd+iz,iob,1)=  &
                              tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)-Nd+iz,iob,1)
  end do
  end do
  end do
end do

!send from kup to kdw

do iob=1,iobnum
  do iz=1,Nd
!$OMP parallel do private(iy,ix)
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_end(3)+iz,iob,1)= &
                              tpsi(ix+mg_sta(1)-1,iy+mg_sta(2)-1,mg_sta(3)+iz-1,iob,1)
  end do
  end do
  end do
end do

end subroutine C_sendrecv_groupob_ngp

!==================================================================================================

end module sendrecv_groupob_ngp_sub
