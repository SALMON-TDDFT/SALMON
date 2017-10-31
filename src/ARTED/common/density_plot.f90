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
subroutine write_density(it,action)
  use Global_Variables
  use salmon_global, only: format3d
  implicit none
  integer :: it
  character(2) :: action
!  character(*) :: position_option

  if(action=='gs') then
    select case(format3d)
    case ('cube')
      write(file_dns_gs,'(2A,"_dns_gs.cube")') trim(directory),trim(SYSname)
      open(502,file=file_dns_gs,position = position_option)
      call write_density_cube(502, .false.)
      close(502)
    case ('vtk')
      write(file_dns_gs,'(2A,"_dns_gs.vtk")') trim(directory),trim(SYSname)
      open(502,file=file_dns_gs,position = position_option)
      call write_density_vtk(502, .false.)
      close(502)
    end select
  endif

  if(action=='rt') then
     select case(format3d)
        case ('cube')
          write(file_dns_rt,200) trim(directory),trim(SYSname),"_dns_rt_", it,".cube"
          write(file_dns_dlt,200)trim(directory),trim(SYSname),"_dns_dlt_",it,".cube"
          open(501,file=file_dns_rt,position = position_option)
          call write_density_cube(501, .false.)
          close(501)
          if(use_ehrenfest_md=='y') then
             call k_shift_wf(Rion_update_rt,Nscf,zu_t,it,"get_dns_gs")
          endif
          open(501,file=file_dns_dlt,position = position_option)
          call write_density_cube(501, .true.)
          close(501)
        case ('vtk')          
          write(file_dns_rt,200) trim(directory),trim(SYSname),"_dns_rt_", it,".vtk"
          write(file_dns_dlt,200)trim(directory),trim(SYSname),"_dns_dlt_",it,".vtk"
          open(501,file=file_dns_rt,position = position_option)
          call write_density_vtk(501, .false.)
          close(501)
          open(501,file=file_dns_dlt,position = position_option)
          call write_density_vtk(501, .true.)
          close(501)
     end select
200  format(3A,I6.6,A)

  endif

end subroutine write_density


subroutine write_density_cube(fh, write_difference)
  use Global_Variables, only: NLx,NLy,NLz,Hx,Hy,Hz,NI,Kion,Rion,Zatom,Lxyz,rho,rho_gs,rho_gs_t,use_ehrenfest_md
  implicit none
  integer, intent(in) :: fh
  logical, intent(in) :: write_difference
  
  integer :: i, ix, iy, iz
  real(8) :: r
    
  write(fh, '(A)') "# SALMON"
  write(fh, '(A)') "# COMMENT"
  write(fh, '(I5,3(F12.6))') NI, 0.00, 0.00, 0.00
  write(fh, '(I5,3(F12.6))') NLx, Hx, 0.00, 0.00
  write(fh, '(I5,3(F12.6))') NLy, 0.00, Hy, 0.00
  write(fh, '(I5,3(F12.6))') NLz, 0.00, 0.00, Hz
  
  do i=1, NI
    write(fh, '(I5,4(F12.6))') Zatom(Kion(i)), 0.00, Rion(1,i), Rion(2,i), Rion(3,i) 
  end do
  
  ! Gaussian .cube file (x-slowest index, z-fastest index)
  i=1
  do ix=0, NLx-1
  do iy=0, NLy-1
  do iz=0, NLz-1
     if (write_difference) then
        if(use_ehrenfest_md=='y') then
           r = rho(Lxyz(ix,iy,iz)) - rho_gs_t(Lxyz(ix,iy,iz))
        else
           r = rho(Lxyz(ix,iy,iz)) - rho_gs(Lxyz(ix,iy,iz))
        endif
     else
        r = rho(Lxyz(ix,iy,iz))
     end if
     if(mod(i,6)==0) then
        write(fh,10) r
     else
        write(fh,10,advance='no') r
     endif
     i=i+1
  end do
  end do
  end do

10 format(ES12.4)
  return
end subroutine write_density_cube


subroutine write_density_vtk(fh, write_difference)
  use Global_Variables, only: NLx, NLy, NLz, Hx, Hy, Hz, Lxyz, rho, rho_gs 
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
        write(fh,10) rho(Lxyz(ix,iy,iz)) - rho_gs(Lxyz(ix,iy,iz))
     else
        write(fh,10) rho(Lxyz(ix,iy,iz))
     end if
  end do
  end do
  end do
10 format(ES12.5)
  return
end subroutine write_density_vtk
  

!! write_density_avs is not inplemented yet!!
!subroutine write_density_avs(fh, write_difference)
!  implicit none
!  integer, intent(in) :: fh
!  logical, intent(in) :: write_difference
!  
!  !todo: please create exporter for "avs express"
!end subroutine write_density_avs




  
  
  
