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
use salmon_parallel, only: nproc_id_global, nproc_group_grid
use salmon_communication, only: comm_is_root, comm_summation
use inputoutput, only: out_dos_start, out_dos_end, out_dos_method, &
                       out_dos_smearing, iout_dos_nenergy, out_dos_fshift, uenergy_from_au
use scf_data
use allocate_psl_sub
use new_world_sub
implicit none
!integer :: iob,iobmax,iob_allob
integer :: iob,iobmax,iob_allob,iik
real(8) :: dos_l_tmp(1:iout_dos_nenergy)
real(8) :: dos_l(1:iout_dos_nenergy)
real(8) :: fk,ww,dw
integer :: iw
real(8) :: ene_homo,ene_lumo,ene_min,ene_max,efermi,eshift

call calc_pmax(iobmax)
ene_min = minval(esp(:,:))
ene_max = maxval(esp(:,:))
ene_homo = esp(nelec/2,1)
ene_lumo = esp(nelec/2+1,1)
if(out_dos_fshift=='y'.and.nstate>nelec/2) then 
  efermi = (ene_homo+ene_lumo)*0.5d0 
  eshift = efermi 
else 
  eshift = 0.d0 
endif 
out_dos_start = max(out_dos_start,ene_min-0.25d0*(ene_max-ene_min))
out_dos_end = min(out_dos_end,ene_max+0.25d0*(ene_max-ene_min))
dw=(out_dos_end-out_dos_start)/dble(iout_dos_nenergy-1) 

dos_l_tmp=0.d0
 
do iik=k_sta,k_end
do iob=1,iobmax
  call calc_allob(iob,iob_allob)
  select case (out_dos_method)
  case('lorentzian') 
    fk=2.d0*out_dos_smearing/pi
    do iw=1,iout_dos_nenergy 
      ww=out_dos_start+dble(iw-1)*dw+eshift-esp(iob_allob,iik)  
      dos_l_tmp(iw)=dos_l_tmp(iw)+wtk(iik)*fk/(ww**2+out_dos_smearing**2) 
    end do 
  case('gaussian')
    fk=2.d0/(sqrt(2.d0*pi)*out_dos_smearing)
    do iw=1,iout_dos_nenergy 
      ww=out_dos_start+dble(iw-1)*dw+eshift-esp(iob_allob,iik)  
      dos_l_tmp(iw)=dos_l_tmp(iw)+wtk(iik)*fk*exp(-(0.5d0/out_dos_smearing**2)*ww**2) 
    end do
  end select
end do
end do
call comm_summation(dos_l_tmp,dos_l,iout_dos_nenergy,nproc_group_grid) 

if(comm_is_root(nproc_id_global))then
  open(101,file="dos.data")
  write(101,'("# Density of States")') 
  select case(unit_energy)
  case('au','a.u.')
    write(101,'("# Energy[a.u.] DOS[a.u.]")') 
  case('ev','eV')
    write(101,'("# Energy[eV]  DOS[1/eV]")') 
  end select
  write(101,'("#-----------------------")') 
  do iw=1,iout_dos_nenergy 
    ww=out_dos_start+dble(iw-1)*dw+eshift
    write(101,'(F16.8,99(1X,E23.15E3))') ww*uenergy_from_au, dos_l(iw)/uenergy_from_au
  end do
  close(101)
end if

end subroutine calc_dos
