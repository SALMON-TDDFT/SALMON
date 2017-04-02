!
!  Copyright 2016 ARTED developers
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
!This file is "Hartree.f90"
!This file contain one subroutine.
!SUBROUTINE Hartree
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130
Subroutine Hartree
  use Global_Variables
  use timer
  implicit none
  integer :: i,ix,iy,iz,n,nx,ny,nz
  real(8) :: G2

  call timer_begin(LOG_HARTREE)

!$omp parallel 

!$omp do private(i)
  do i=1,NL
    rho_3D(Lx(i),Ly(i),Lz(i))=rho(i)
  end do
!$omp end do

!$omp do private(ix,iy,nz) collapse(3)
  do nz = -NLz/2,NLz/2-1
  do iy = 0,NLy-1
  do ix = 0,NLx-1
    f1(ix,iy,nz)=sum(eGzc(nz,:)*rho_3D(ix,iy,:))
  end do
  end do
  end do
!$omp end do

!$omp do private(ix,ny,nz) collapse(3)
  do nz = -NLz/2,NLz/2-1
  do ny = -NLy/2,NLy/2-1
  do ix = 0,NLx-1
    f2(ix,ny,nz)=sum(eGyc(ny,:)*f1(ix,:,nz))
  end do
  end do
  end do
!$omp end do

!$omp do private(nx,ny,nz) collapse(3)
  do nz = -NLz/2,NLz/2-1
  do ny = -NLy/2,NLy/2-1
  do nx = -NLx/2,NLx/2-1
    rhoe_G_3D(nx,ny,nz)=sum(eGxc(nx,:)*f2(:,ny,nz))/dble(NL)
  end do
  end do
  end do
!$omp end do

!$omp do private(nx,ny,nz,G2) collapse(3)
  do nz = -NLz/2,NLz/2-1
  do ny = -NLy/2,NLy/2-1
  do nx = -NLx/2,NLx/2-1
    rhoe_G_temp(nxyz(nx,ny,nz))=rhoe_G_3D(nx,ny,nz)
    G2=Gx(nxyz(nx,ny,nz))**2+Gy(nxyz(nx,ny,nz))**2+Gz(nxyz(nx,ny,nz))**2
    rhoe_G_3D(nx,ny,nz)=4*Pi/G2*rhoe_G_3D(nx,ny,nz)
  end do
  end do
  end do
!$omp end do

  rhoe_G_3D(0,0,0)=0.d0

!$omp do private(n)
  do n=NG_s,NG_e
    rhoe_G(n)=rhoe_G_temp(n)
  end do
!$omp end do

!$omp do private(nx,ny,iz) collapse(3)
  do iz = 0,NLz-1
  do ny = -NLy/2,NLy/2-1
  do nx = -NLx/2,NLx/2-1
    f3(nx,ny,iz)=sum(eGz(:,iz)*rhoe_G_3D(nx,ny,:))
  end do
  end do
  end do
!$omp end do

!$omp do private(nx,iy,iz) collapse(3)
  do iz = 0,NLz-1
  do iy = 0,NLy-1
  do nx = -NLx/2,NLx/2-1
    f4(nx,iy,iz)=sum(eGy(:,iy)*f3(nx,:,iz))
  end do
  end do
  end do
!$omp end do

!$omp do private(ix,iy,iz) collapse(3)
  do iz = 0,NLz-1
  do iy = 0,NLy-1
  do ix = 0,NLx-1
    Vh_3D(ix,iy,iz)=sum(eGx(:,ix)*f4(:,iy,iz))
  end do
  end do
  end do
!$omp end do

!$omp do private(i)
  do i=1,NL
    Vh(i)=Vh_3D(Lx(i),Ly(i),Lz(i))
  end do
!$omp end do

!$omp end parallel

  call timer_end(LOG_HARTREE)

  return
End Subroutine Hartree
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------130

