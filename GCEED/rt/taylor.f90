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
subroutine taylor(tzpsi_in,tzpsi_out,htpsi)
use scf_data
use allocate_mat_sub
use deallocate_mat_sub
implicit none

integer :: nn,ix,iy,iz
complex(8) :: tzpsi_in(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                   mg_sta(2)-Nd:mg_end(2)+Nd,    &
                   mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
complex(8) :: htpsi(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                    mg_sta(2)-Nd:mg_end(2)+Nd,    &
                    mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)
complex(8) :: tzpsi_out(mg_sta(1)-Nd:mg_end(1)+Nd+1,    &
                   mg_sta(2)-Nd:mg_end(2)+Nd,    &
                   mg_sta(3)-Nd:mg_end(3)+Nd,1:iobnum,1)

!call MPI_BARRIER(nproc_group_orbital,ierr)

iwk_size=2
call make_iwksta_iwkend

if(ilsda==0)then
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhobox(ix,iy,iz) = 0.d0
  end do
  end do
  end do
else if(ilsda==1)then
!$OMP parallel do
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    rhobox_s(ix,iy,iz,1) = 0.d0
    rhobox_s(ix,iy,iz,2) = 0.d0
  end do
  end do
  end do
end if

if(ihpsieff==1)then
  if(ilsda==0)then
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
    end do
    end do
    end do
  else if(ilsda==1)then 
!$OMP parallel do
    do iz=mg_sta(3),mg_end(3)
    do iy=mg_sta(2),mg_end(2)
    do ix=mg_sta(1),mg_end(1)
      Vlocal2(ix,iy,iz,1)=Vlocal(ix,iy,iz,1)+Vbox(ix,iy,iz)
      Vlocal2(ix,iy,iz,2)=Vlocal(ix,iy,iz,2)+Vbox(ix,iy,iz)
    end do
    end do
    end do
  end if
end if

!  call MPI_BARRIER(nproc_group_orbital,ierr)

  do nn=1,N_hamil
    if(ihpsieff==1)then
      if(mod(nn,2)==1)then
        call hpsi_groupob(tzpsi_in,htpsi,tzpsi_out,Vlocal2,nn,1)
      else
        call hpsi_groupob(htpsi,tzpsi_in,tzpsi_out,Vlocal2,nn,1)
      end if
    else
      if(mod(nn,2)==1)then
        call hpsi_groupob(tzpsi_in,htpsi,tzpsi_out,Vlocal,nn,1)
      else
        call hpsi_groupob(htpsi,tzpsi_in,tzpsi_out,Vlocal,nn,1)
      end if
    end if
  end do

!call MPI_BARRIER(nproc_group_orbital,ierr)

end subroutine taylor

