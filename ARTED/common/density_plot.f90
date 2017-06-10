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

subroutine write_density_cube(fh, write_difference)
  use Global_Variables, only: NLx, NLy, NLz, aLx, aLy, aLz, NI, Kion, Rion, Zatom, Lxyz, Rho, Rho_gs 
  implicit none
  integer, intent(in) :: fh
  logical, intent(in) :: write_difference
  
  integer :: i, ix, iy, iz
    
  write(fh, '(A)') "# SALMON"
  write(fh, '(A)') "# COMMENT"
  write(fh, '(I2,3(1X,F5.2))') NI, 0.00, 0.00, 0.00
  write(fh, '(I2,3(1X,F5.2))') NLx, aLx, 0.00, 0.00
  write(fh, '(I2,3(1X,F5.2))') NLy, 0.00, aLy, 0.00
  write(fh, '(I2,3(1X,F5.2))') NLz, 0.00, 0.00, aLz
  
  do i=1, NI
    write(fh, '(I2,4(1X,F5.2))') Zatom(Kion(i)), 0.00, Rion(1,i)*aLx, Rion(2,i)*aLy, Rion(3,i)*aLz 
  end do
  
  ! Gaussian .cube file (x-slowest index, z-fastest index)
  do ix=0, NLx-1
    do iy=0, NLy-1
      do iz=0, NLz-1
        if (write_difference) then
          write(fh, '(ES12.5)') Rho(Lxyz(ix, iy, iz)) - Rho_gs(Lxyz(ix, iy, iz))
        else
          write(fh, '(ES12.5)') Rho(Lxyz(ix, iy, iz))
        end if
      end do
    end do
  end do
end subroutine write_density_cube



subroutine write_density_vtk(fh, write_difference)
  use Global_Variables, only: NLx, NLy, NLz, aLx, aLy, aLz, Hx, Hy, Hz, Lxyz, Rho, Rho_gs 
  implicit none
  integer, intent(in) :: fh
  logical, intent(in) :: write_difference
  integer :: ix, iy, iz
  
  write(fh, '(A)') "# vtk DataFile Version 3.0"
  write(fh, '(A)') "vtk output"
  write(fh, '(A)') "ASCII"
  write(fh, '(A)') "DATASET STRUCTURED_POINTS"
  write(fh, '(A,3(1X,I2))') "DIMENSIONS", NLx, NLy, NLz
  write(fh, '(A,3(1X,F3.1))') "ORIGIN", 0.0, 0.0, 0.0
  write(fh, '(A,3(1X,F7.3))') "SPACING", Hx, Hy, Hz
  write(fh, '(A,1X,I6)') "POINT_DATA", NLx * NLy * NLz
  write(fh, '(A)') "SCALARS scalars float"
  write(fh, '(A)') "LOOKUP_TABLE default"
  
  ! VTK file (x-fastest index, z-slowest index)
  do iz=0, NLz-1
    do iy=0, NLy-1
      do ix=0, NLx-1
        if (write_difference) then
          write(fh, '(ES12.5)') Rho(Lxyz(ix, iy, iz)) - Rho_gs(Lxyz(ix, iy, iz))
        else
          write(fh, '(ES12.5)') Rho(Lxyz(ix, iy, iz))
        end if
      end do
    end do
  end do
  return
end subroutine write_density_vtk
  


subroutine write_density_avs(fh, write_difference)
  implicit none
  integer, intent(in) :: fh
  logical, intent(in) :: write_difference
  
  !todo: please create exporter for "avs express"
end subroutine write_density_avs




  
  
  
