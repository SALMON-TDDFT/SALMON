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

subroutine writeDensity(cOutFileName,iwdenoption2,iwdenstep2,denplane2,idensum2,posplane2,iplType)
use scf_data
implicit none
character(30),intent(in):: cOutFileName
character(30):: cOutFileName2
integer::ix,iy,iz
integer::jj
integer::jsta,jend
real(8),allocatable::rRho(:,:,:) 
real(8),allocatable::rRho2(:,:,:)
integer::iplType
integer       :: iwdenoption2
integer       :: iwdenstep2
character(2)  :: denplane2 
integer       :: idensum2   
real(8)       :: posplane2
character(8)  :: fileNumber_data
integer::icount

allocate( rRho(lg_sta(1):lg_end(1), lg_sta(2):lg_end(2), lg_sta(3):lg_end(3)) )
allocate( rRho2(lg_sta(1):lg_end(1), lg_sta(2):lg_end(2), lg_sta(3):lg_end(3)) )
if(iplType == 0) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),   &
       ng_sta(2):ng_end(2),   &
       ng_sta(3):ng_end(3))   &
    = rho(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3))
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rRho2 = rRho2/a_B**3
elseif(iplType == 1) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),   &
       ng_sta(2):ng_end(2),   &
       ng_sta(3):ng_end(3))   &
    = rho(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3))   &
    - rho0(ng_sta(1):ng_end(1),   &
           ng_sta(2):ng_end(2),   &
           ng_sta(3):ng_end(3))
       
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rRho2 = rRho2/a_B**3
elseif(iplType == 2) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),   &
       ng_sta(2):ng_end(2),   &
       ng_sta(3):ng_end(3))   &
    = elf(ng_sta(1):ng_end(1),   &
          ng_sta(2):ng_end(2),   &
          ng_sta(3):ng_end(3))
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
elseif(iplType == 10) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3))   &
    = dble(Ex_static(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rRho2 = rRho2*5.14223d1  !  unit is Volt/Angstrom
elseif(iplType == 11) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3))   &
    = dble(Ey_static(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rRho2 = rRho2*5.14223d1  !  unit is Volt/Angstrom
elseif(iplType == 12) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3))   &
    = dble(Ez_static(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rRho2 = rRho2*5.14223d1  !  unit is Volt/Angstrom
elseif(iplType == 13) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3))   &
    = dble(Eeff_dip(1,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rRho2 = rRho2*5.14223d1  !  unit is Volt/Angstrom
elseif(iplType == 14) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3))   &
    = dble(Eeff_dip(2,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rRho2 = rRho2*5.14223d1  !  unit is Volt/Angstrom
elseif(iplType == 15) then
  rRho=0.d0
  rRho(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3))   &
    = dble(Eeff_dip(3,ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
  call MPI_Allreduce(rRho,rRho2, &
                   lg_num(1)*lg_num(2)*lg_num(3), &
                   MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  rRho2 = rRho2*5.14223d1  !  unit is Volt/Angstrom
endif

if(numfile_movie>=2.and.iwdenoption2==1)then
  if(myrank<numfile_movie)then
    write(fileNumber_data, '(i8)') myrank
    cOutFileName2 = trim(cOutFileName)//"."//adjustl(fileNumber_data)
    open(1,file=cOutFileName2)
    jsta=myrank*(lg_num(1)*lg_num(2)*lg_num(3))/numfile_movie+1
    jend=(myrank+1)*(lg_num(1)*lg_num(2)*lg_num(3))/numfile_movie
    do jj=jsta,jend
      if(abs(rRho2(icoo1d(1,jj),icoo1d(2,jj),icoo1d(3,jj)))>=1.0d-6)then
        write(1,'(e20.8)') rRho2(icoo1d(1,jj),icoo1d(2,jj),icoo1d(3,jj))
      else
        write(1,'(a1)') "0"
      end if
    enddo
    close(1)
  end if
else
  if(myrank==0)then
    cOutFileName2=trim(adjustl(cOutFileName))//".dx"
    open(1,file=cOutFileName2)
    if(iwdenoption2==1)then
      call output_dx_header(1)
      do ix=lg_sta(1),lg_end(1),1
      do iy=lg_sta(2),lg_end(2),1
      do iz=lg_sta(3),lg_end(3),1
        if(abs(rRho2(ix,iy,iz))>=1.0d-10)then
          write(1,'(e20.8)') rRho2(ix,iy,iz)
        else
          write(1,'(a1)') "0"
        end if
      enddo
      enddo
      enddo
    end if
    close(1)
  endif
endif

deallocate ( rRho, rRho2 )

end subroutine

SUBROUTINE output_dx_header(fp)

use scf_data
use allocate_mat_sub
implicit none

integer, intent(IN) :: fp

integer :: Ngridx, Ngridy, Ngridz
real(8) :: originx, originy, originz, dx, dy, dz

Ngridx = mg_end(1)-mg_sta(1)+1
Ngridy = mg_end(2)-mg_sta(2)+1
Ngridz = mg_end(3)-mg_sta(3)+1

originx = vecR(1, mg_sta(1), mg_sta(2), mg_sta(3))*Hgs(1)*a_B
originy = vecR(2, mg_sta(1), mg_sta(2), mg_sta(3))*Hgs(2)*a_B
originz = vecR(3, mg_sta(1), mg_sta(2), mg_sta(3))*Hgs(3)*a_B

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

END SUBROUTINE output_dx_header
