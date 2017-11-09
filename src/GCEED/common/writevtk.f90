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
subroutine writevtk(fp,suffix,tmatbox_l)
  use inputoutput, only: au_length_aa
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use scf_data
  implicit none
  integer, intent(IN) :: fp
  character(30),intent(in):: suffix
  real(8), intent(IN) :: tmatbox_l(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),  &
                                   lg_sta(3):lg_end(3))
  character(30):: filename
  integer :: ix,iy,iz

  if(comm_is_root(nproc_id_global))then
    filename=trim(suffix)//".vtk"
    open(fp,file=filename)

    write(fp, '(A)') "# vtk DataFile Version 3.0"
    write(fp, '(A)') "vtk output"
    write(fp, '(A)') "ASCII"
    write(fp, '(A)') "DATASET STRUCTURED_POINTS"
    write(fp, '(A,3(1X,I2))') "DIMENSIONS", lg_num(1), lg_num(2), lg_num(3)
    write(fp, '(A,3(1X,F10.5))') "ORIGIN",gridcoo(lg_sta(1),1)*au_length_aa,  &
                                          gridcoo(lg_sta(2),2)*au_length_aa,  &
                                          gridcoo(lg_sta(3),3)*au_length_aa
    write(fp, '(A,3(1X,F10.5))') "SPACING", Hgs(1)*au_length_aa,  &
                                            Hgs(2)*au_length_aa,  &
                                            Hgs(3)*au_length_aa
    write(fp, '(A,1X,I6)') "POINT_DATA", lg_num(1) * lg_num(2) * lg_num(3)
    write(fp, '(A)') "SCALARS scalars float"
    write(fp, '(A)') "LOOKUP_TABLE default"


    do iz=lg_sta(3),lg_end(3)
    do iy=lg_sta(2),lg_end(2)
    do ix=lg_sta(1),lg_end(1)
      write(fp, '(ES12.5)') tmatbox_l(ix,iy,iz)
    end do
    end do
    end do
  
    close(fp)
  end if

end subroutine writevtk

