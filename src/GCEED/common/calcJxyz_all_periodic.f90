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
subroutine calcJxyz_all_periodic
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation
  use scf_data
  use read_pslfile_sub
  use allocate_psl_sub
  implicit none
  
  integer :: ix,iy,iz
  integer :: iatom
  
  complex(8),parameter :: zI=(0.d0,1.d0)
  
  integer :: i,j
  integer :: mmx(mg_num(1)*mg_num(2)*mg_num(3))
  integer :: mmy(mg_num(1)*mg_num(2)*mg_num(3))
  integer :: mmz(mg_num(1)*mg_num(2)*mg_num(3))
  integer :: lx(lg_num(1)*lg_num(2)*lg_num(3))
  integer :: ly(lg_num(1)*lg_num(2)*lg_num(3))
  integer :: lz(lg_num(1)*lg_num(2)*lg_num(3))
  real(8) :: alx,aly,alz
  real(8) :: hx,hy,hz

! nonlocal potential
  if(comm_is_root(nproc_id_global))then
    write(*,*) ''
    write(*,*) '============nonlocal grid data=============='
  endif

  hx=Hgs(1) 
  hy=Hgs(2) 
  hz=Hgs(3)
  alx=Hgs(1)*dble(lg_num(1))
  aly=Hgs(2)*dble(lg_num(2))
  alz=Hgs(3)*dble(lg_num(3))

  do iz=1,mg_num(3)
  do iy=1,mg_num(2)
  do ix=1,mg_num(1)
    i=(iz-1)*mg_num(1)*mg_num(2)+(iy-1)*mg_num(1)+ix
    mmx(i)=ix+mg_sta(1)-1
    mmy(i)=iy+mg_sta(2)-1
    mmz(i)=iz+mg_sta(3)-1
  end do
  end do
  end do
 
  do iz=1,lg_num(3)
  do iy=1,lg_num(2)
  do ix=1,lg_num(1)
    i=(iz-1)*lg_num(1)*lg_num(2)+(iy-1)*lg_num(1)+ix
    lx(i)=ix
    ly(i)=iy
    lz(i)=iz
  end do
  end do
  end do
 
  call calc_mps(pp,ppg,alx,aly,alz,lx,ly,lz,lg_num(1)*lg_num(2)*lg_num(3),   &
                                   mmx,mmy,mmz,mg_num(1)*mg_num(2)*mg_num(3),   &
                                   hx,hy,hz)
  call calc_mps(pp,ppg_all,alx,aly,alz,lx,ly,lz,lg_num(1)*lg_num(2)*lg_num(3),   &
                                       lx,ly,lz,lg_num(1)*lg_num(2)*lg_num(3),   &
                                       hx,hy,hz)
  Mps(1:MI)=ppg%mps(1:MI) 
  Mps_all(1:MI)=ppg_all%mps(1:MI) 

  call init_jxyz(ppg) 
  call init_jxyz(ppg_all)
 
  call calc_jxyz(pp,ppg,alx,aly,alz,lx,ly,lz,lg_num(1)*lg_num(2)*lg_num(3),   &
                                    mmx,mmy,mmz,mg_num(1)*mg_num(2)*mg_num(3),   &
                                    hx,hy,hz)
  call calc_jxyz(pp,ppg_all,alx,aly,alz,lx,ly,lz,lg_num(1)*lg_num(2)*lg_num(3),   &
                                    lx,ly,lz,lg_num(1)*lg_num(2)*lg_num(3),   &
                                    hx,hy,hz)
  Jxyz=0
  Jxyz_all=0
  Jxxyyzz=0
  Jxxyyzz_all=0

  do iatom=1,MI
    do j=1,Mps(iatom)
      Jxyz(1,j,iatom)=ppg%jxyz(1,j,iatom)
      Jxyz(2,j,iatom)=ppg%jxyz(2,j,iatom)
      Jxyz(3,j,iatom)=ppg%jxyz(3,j,iatom)
      Jxxyyzz(1,j,iatom)=ppg%jxx(j,iatom)
      Jxxyyzz(2,j,iatom)=ppg%jyy(j,iatom)
      Jxxyyzz(3,j,iatom)=ppg%jzz(j,iatom)
    end do
    do j=1,Mps_all(iatom)
      Jxyz_all(1,j,iatom)=ppg_all%jxyz(1,j,iatom)
      Jxyz_all(2,j,iatom)=ppg_all%jxyz(2,j,iatom)
      Jxyz_all(3,j,iatom)=ppg_all%jxyz(3,j,iatom)
      Jxxyyzz_all(1,j,iatom)=ppg_all%jxx(j,iatom)
      Jxxyyzz_all(2,j,iatom)=ppg_all%jyy(j,iatom)
      Jxxyyzz_all(3,j,iatom)=ppg_all%jzz(j,iatom)
    end do
  end do
  
  if(iSCFRT==1)then
    do iatom=1,MI
      if(comm_is_root(nproc_id_global))then
        write(*,*) "Mps =", Mps_all(iatom)
      end if
    end do
  end if
  
  return

end subroutine calcJxyz_all_periodic
