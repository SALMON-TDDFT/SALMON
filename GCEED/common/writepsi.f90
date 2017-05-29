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
!======================================================================
!======================================================================
subroutine writepsi
use scf_data
use allocate_mat_sub
implicit none
integer :: iob,ix,iy,iz
integer :: p0,iob_myob,icheck_corrkob
character(50) :: filePsi
character(10) :: filenum

if(iSCFRT==1)then
  do p0=1,itotMST
    write(filenum, '(i3)') p0
    filePsi = "psi"//trim(adjustl(filenum))//".cube"
    if(myrank==0) then
      open(103,file=filePsi)
      call output_dx_header_psi(103)
    end if
    call conv_p0(p0,iob)
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,icheck_corrkob)
!OMP parallel do
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do
    if(icheck_corrkob==1)then
!OMP parallel do
      do iz=mg_sta(3),lg_end(3)
      do iy=mg_sta(2),lg_end(2)
      do ix=mg_sta(1),lg_end(1)
        matbox_l(ix,iy,iz)=psi(ix,iy,iz,iob_myob,1)
      end do
      end do
      end do
    end if
    call MPI_Allreduce(matbox_l,matbox_l2,  &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myrank==0)then
      do ix=lg_sta(1),lg_end(1)
      do iy=lg_sta(2),lg_end(2)
      do iz=lg_sta(3),lg_end(3)
        if(mod(iz+1-lg_sta(3),6)==0)then
          write(103,'(e13.5)', advance="yes") matbox_l2(ix,iy,iz)
        else
          write(103,'(e13.5)', advance="no") matbox_l2(ix,iy,iz)
        endif 
      end do 
      write(103,*) 
      end do 
      end do 
      close(103)
    end if
  end do
else if(iSCFRT==2)then
  do p0=1,itotMST
    write(filenum, '(i3)') p0
    filePsi = "psi"//trim(adjustl(filenum))//".cube"
    if(myrank==0)then
      open(103,file=filePsi)
      call output_dx_header_psi(103)
    end if
    call conv_p0(p0,iob)
    call calc_myob(iob,iob_myob)
    call check_corrkob(iob,icheck_corrkob)
!OMP parallel do
    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      cmatbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do
    if(icheck_corrkob==1)then
!OMP parallel do
      do iz=mg_sta(3),lg_end(3)
      do iy=mg_sta(2),lg_end(2)
      do ix=mg_sta(1),lg_end(1)
        cmatbox_l(ix,iy,iz)=zpsi(ix,iy,iz,iob_myob,1)
      end do
      end do
      end do
    end if
    call MPI_Allreduce(cmatbox_l,cmatbox_l2,  &
&             lg_num(1)*lg_num(2)*lg_num(3), &
&             MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(myrank==0)then
      do ix=lg_sta(1),lg_end(1)
      do iy=lg_sta(2),lg_end(2)
      do iz=lg_sta(3),lg_end(3)
        if(mod(iz+1-lg_sta(3),6)==0)then
          write(103,'(e13.5)', advance="yes") abs(cmatbox_l2(ix,iy,iz))
        else
          write(103,'(e13.5)', advance="no") abs(cmatbox_l2(ix,iy,iz))
        endif 
      end do 
      write(103,*) 
      end do 
      end do 
      close(103)
    end if
  end do
end if

end subroutine writepsi

SUBROUTINE output_dx_header_psi(fp)

use scf_data
use allocate_mat_sub
implicit none

integer, intent(IN) :: fp

integer :: Ngridx, Ngridy, Ngridz
real(8) :: originx, originy, originz, dx, dy, dz

Ngridx = lg_end(1)-lg_sta(1)+1
Ngridy = lg_end(2)-lg_sta(2)+1
Ngridz = lg_end(3)-lg_sta(3)+1

if (imesh_oddeven == 1) then
    originx = dble(lg_sta(1))*Hgs(1)*a_B
    originy = dble(lg_sta(2))*Hgs(2)*a_B
    originz = dble(lg_sta(3))*Hgs(3)*a_B
else if (imesh_oddeven == 2) then
    originx = dble(lg_sta(1))*Hgs(1)*a_B
    originy = dble(lg_sta(2))*Hgs(2)*a_B
    originz = dble(lg_sta(3))*Hgs(3)*a_B
end if

dx = Hgs(1)*a_B
dy = Hgs(2)*a_B
dz = Hgs(3)*a_B

write(fp, '(a)', advance = "no") "object 1 class gridpositions counts"
write(fp, '(3i5)', advance = "yes") Ngridx, Ngridy, Ngridz
write(fp, '(a7)', advance = "no") " origin"
write(fp, '(3f12.6)', advance = "yes") originx, originy, originz
write(fp, '(a7)', advance = "no") " delta "
write(fp, '(3f12.6)', advance = "yes") dx, 0.d0, 0.d0
write(fp, '(a7)', advance = "no") " delta "
write(fp, '(3f12.6)', advance = "yes") 0.d0, dy, 0.d0
write(fp, '(a7)', advance = "no") " delta "
write(fp, '(3f12.6)', advance = "yes") 0.d0, 0.d0, dz
write(fp, '(a)', advance = "no") "object 2 class gridconections counts"
write(fp, '(3i5)', advance = "yes") Ngridx, Ngridy, Ngridz
write(fp, '(a)', advance = "no") "object 3 class array type float rank 0 items"
write(fp, '(i8)', advance = "no") Ngridx*Ngridy*Ngridz
write(fp, '(a20)', advance = "yes") "data follows"

END SUBROUTINE output_dx_header_psi
