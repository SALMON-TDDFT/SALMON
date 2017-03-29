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

subroutine copy_density
use scf_data
implicit none
integer :: iiter
integer :: is

if(Miter==1)then
  rho_in(:,:,:,num_rho_stock+1)=rho(:,:,:)
  if(ilsda==1)then
    rho_s_in(:,:,:,1:2,num_rho_stock+1)=rho_s(:,:,:,1:2)
  end if
end if

do iiter=1,num_rho_stock
  rho_in(:,:,:,iiter)=rho_in(:,:,:,iiter+1)
end do
do iiter=1,num_rho_stock-1
  rho_out(:,:,:,iiter)=rho_out(:,:,:,iiter+1)
end do

if(ilsda==1)then
  do iiter=1,num_rho_stock
    do is=1,2
      rho_s_in(:,:,:,is,iiter)=rho_s_in(:,:,:,is,iiter+1)
    end do
  end do
  do iiter=1,num_rho_stock-1
    do is=1,2
      rho_s_out(:,:,:,is,iiter)=rho_s_out(:,:,:,is,iiter+1)
    end do
  end do
end if

end subroutine copy_density
