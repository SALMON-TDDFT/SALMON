!
!  Copyright 2018 SALMON developers
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
subroutine buffer_broyden_ns(iter)
  use salmon_parallel, only: nproc_group_global
  use broyden_sub
  use scf_data, only: ng_sta,ng_end,ng_num,num_rho_stock,rho,rho_in,rho_out,  &
                      ilsda,rho_s,rho_s_in,rho_s_out,MST,ifMST
  implicit none
  integer,intent(in) :: iter
  integer :: ix,iy,iz,is
  real(8) :: vecr(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3))

  if(ilsda==0)then

    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      vecr(ix,iy,iz)=rho(ix,iy,iz)
    end do
    end do
    end do

    call broyden(vecr,rho_in,rho_out,ng_num(1)*ng_num(2)*ng_num(3),  &
                 iter,num_rho_stock,num_rho_stock,nproc_group_global)
  
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      rho(ix,iy,iz)= vecr(ix,iy,iz)
    end do
    end do
    end do

  else if(ilsda==1)then
    
    do is=1,2
      if(ifMST(is)>=1.and.MST(is)>=1)then
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          vecr(ix,iy,iz)=rho_s(ix,iy,iz,is)
        end do
        end do
        end do
  
        call broyden(vecr,rho_s_in(ng_sta(1):,ng_sta(2):,ng_sta(3):,1:,is),  &
                     rho_s_out(ng_sta(1):,ng_sta(2):,ng_sta(3):,1:,is),  &
                     ng_num(1)*ng_num(2)*ng_num(3),  &
                     iter,num_rho_stock,num_rho_stock,nproc_group_global)
  
        do iz=ng_sta(3),ng_end(3)
        do iy=ng_sta(2),ng_end(2)
        do ix=ng_sta(1),ng_end(1)
          rho_s(ix,iy,iz,is)= vecr(ix,iy,iz)
        end do
        end do
        end do
      end if
    end do

    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      rho(ix,iy,iz)= rho_s(ix,iy,iz,1)+rho_s(ix,iy,iz,2)
    end do
    end do
    end do
  end if

end subroutine buffer_broyden_ns
