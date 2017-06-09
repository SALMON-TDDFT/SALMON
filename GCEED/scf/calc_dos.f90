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
subroutine calc_dos
use inputoutput
use scf_data
use allocate_psl_sub
use new_world_sub
implicit none
integer :: iob,iobmax,iob_allob,iene
real(8) :: rbox_dos(-300:300)
real(8) :: dos(-300:300)
real(8),parameter :: sigma_gd=0.01d0

call calc_pmax(iobmax)

rbox_dos=0.d0

do iob=1,iobmax
  call calc_allob(iob,iob_allob)
  do iene=-300,300
    rbox_dos(iene)=rbox_dos(iene)  &
              +exp(-(dble(iene)/10d0/au_energy_ev-esp(iob_allob,1))**2/(2.d0*sigma_gd**2))/sqrt(2.d0*Pi*sigma_gd**2)
  end do
end do
call MPI_Allreduce(rbox_dos,dos,601,MPI_DOUBLE_PRECISION,MPI_SUM,newworld_comm_grid,ierr) 

if(myrank==0)then
  open(101,file="dos.data")
  do iene=-300,300
    write(101,'(f10.5,f14.8)') dble(iene)/10.d0/au_energy_ev*uenergy_from_au,&
                               dos(iene)*au_energy_ev/uenergy_from_au
  end do
  close(101)
end if

end subroutine calc_dos
