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
subroutine writecube(fp,suffix,phys_quantity,tmatbox_l)
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use scf_data
  implicit none
  integer, intent(IN) :: fp
  character(30),intent(in):: suffix
  character(30),intent(in):: phys_quantity
  real(8), intent(IN) :: tmatbox_l(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),  &
                                   lg_sta(3):lg_end(3))
  character(30):: filename
  integer :: j,iatom
  integer :: ix,iy,iz
  integer :: ik

  if(comm_is_root(nproc_id_global))then
    filename=trim(suffix)//".cube"
    open(fp,file=filename)
    if(phys_quantity=="psi")then
      write(fp,*) "Molecular Orbital"
    else if(phys_quantity=="dns")then
      write(fp,*) "Electron Density"
    else if(phys_quantity=="elf")then
      write(fp,*) "Electron Localization Function"
    else if(phys_quantity=="dnsdiff")then
      write(fp,*) "Difference of Electron Density"
    else if(phys_quantity=="exsta")then
      write(fp,*) "x Component of Static Electric Field"
    else if(phys_quantity=="eysta")then
      write(fp,*) "y Component of Static Electric Field"
    else if(phys_quantity=="ezsta")then
      write(fp,*) "z Component of Static Electric Field"
    end if
    write(fp,*) "All values here are in a.u."
    write(fp,'(i5,3f12.6)') MI,gridcoo(lg_sta(1),1),gridcoo(lg_sta(2),2),gridcoo(lg_sta(3),3)
    write(fp,'(i5,3f12.6)') lg_num(1),Hgs(1),0.d0,0.d0
    write(fp,'(i5,3f12.6)') lg_num(2),0.d0,Hgs(2),0.d0
    write(fp,'(i5,3f12.6)') lg_num(3),0.d0,0.d0,Hgs(3)
    do iatom=1,MI
      ik=Kion(iatom)
      write(fp,'(i5,4f12.6)') iZatom(ik),dble(iZatom(ik)),(Rion(j,iatom),j=1,3)
    end do

    do ix=lg_sta(1),lg_end(1)
    do iy=lg_sta(2),lg_end(2)
      write(fp,'(6e13.5)', advance="yes") (tmatbox_l(ix,iy,iz),iz=lg_sta(3),lg_end(3))
!    do iz=lg_sta(3),lg_end(3)
!      if(mod(iz+1-lg_sta(3),6)==0)then
!        write(fp,'(e13.5)', advance="yes") abs(tmatbox_l(ix,iy,iz))
!      else
!        write(fp,'(e13.5)', advance="no") abs(tmatbox_l(ix,iy,iz))
!      endif
!    end do
!    write(fp,*)
    end do
    end do
    close(fp)
  end if

end subroutine writecube

