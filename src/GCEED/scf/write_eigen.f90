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
subroutine write_eigen
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root
use inputoutput, only: uenergy_from_au
use scf_data
implicit none
integer :: iob

if(comm_is_root(nproc_id_global))then
  open(101,file=file_eigen)
  write(101,'("# 1 particle energies")') 
  select case(unit_energy)
  case('au','a.u.')
    write(101,'("# Orbital   Energy[a.u.]")') 
  case('ev','eV')
    write(101,'("# Orbital   Energy[eV]")') 
  end select
  write(101,'("#-----------------------")') 
  do iob=1,itotmst
    write(101,'(1x,i5,e26.16e3)') iob, esp(iob,1)*uenergy_from_au
  end do
  close(101)
end if

end subroutine write_eigen
