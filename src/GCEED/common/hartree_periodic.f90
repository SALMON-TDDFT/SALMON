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
!SUBROUTINE Hartree_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------
subroutine Hartree_periodic(trho,tVh)
  use salmon_parallel, only: nproc_group_global, nproc_group_bound
  use salmon_communication, only: comm_summation
  use scf_data
  use new_world_sub
  use allocate_mat_sub
  implicit none
  integer :: ix,iy,iz,kx,ky,kz,kkx,kky,kkz
  real(8) :: trho(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  real(8) :: tVh(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
  real(8) :: Gx,Gy,Gz
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8) :: G2
  real(8) :: bLx,bLy,bLz
  integer :: n

!$OMP parallel do private(iz,iy,ix)
  do iz=lg_sta(3),lg_end(3)
  do iy=lg_sta(2),lg_end(2)
  do ix=lg_sta(1),lg_end(1)
    ff1(ix,iy,iz)=0.d0
  end do
  end do
  end do
!$OMP parallel do
  do n=1,lg_num(1)*lg_num(2)*lg_num(3)
    rhoe_G_tmp(n)=0.d0
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=lg_sta(3),lg_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    trho2z(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    trho2z(ix,iy,iz)=trho(ix,iy,iz)
  end do
  end do
  end do

  call comm_summation(trho2z,trho3z,ng_num(1)*ng_num(2)*lg_num(3),nproc_group_bound(3))
  
!$OMP parallel do private(ix,kx)
  do ix=lg_sta(1),lg_end(1)
    do kx=lg_sta(1),lg_end(1)
      eGx(kx,ix)=exp(zI*(2.d0*Pi*dble((ix-1)*(kx-1))/dble(lg_num(1))))
      eGxc(kx,ix)=conjg(eGx(kx,ix))
    end do
  end do
!$OMP parallel do private(iy,ky)
  do iy=lg_sta(2),lg_end(2)
    do ky=lg_sta(2),lg_end(2)
      eGy(ky,iy)=exp(zI*(2.d0*Pi*dble((iy-1)*(ky-1))/dble(lg_num(2))))
      eGyc(ky,iy)=conjg(eGy(ky,iy))
    end do
  end do
!$OMP parallel do private(iz,kz)
  do iz=lg_sta(3),lg_end(3)
    do kz=lg_sta(3),lg_end(3)
      eGz(kz,iz)=exp(zI*(2.d0*Pi*dble((iz-1)*(kz-1))/dble(lg_num(3))))
      eGzc(kz,iz)=conjg(eGz(kz,iz))
    end do
  end do

!$OMP parallel do private(iz,iy,ix)
  do iz=ng_sta(3),ng_end(3)
  do iy=lg_sta(2),lg_end(2)
  do ix=ng_sta(1),ng_end(1)
    ff1y(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,iy,ix)
  do kz = ng_sta(3),ng_end(3)
  do iy = ng_sta(2),ng_end(2)
  do ix = ng_sta(1),ng_end(1)
    ff1y(ix,iy,kz)=sum(eGzc(kz,:)*trho3z(ix,iy,:))
  end do
  end do
  end do
  call comm_summation(ff1y,ff2y,ng_num(1)*lg_num(2)*ng_num(3),nproc_group_bound(2))

!$OMP parallel do private(iz,iy,ix)
  do iz=ng_sta(3),ng_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=lg_sta(1),lg_end(1)
    ff1x(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,ky,ix)
  do kz = ng_sta(3),ng_end(3)
  do ky = ng_sta(2),ng_end(2)
  do ix = ng_sta(1),ng_end(1)
    ff1x(ix,ky,kz)=sum(eGyc(ky,:)*ff2y(ix,:,kz))
  end do
  end do
  end do

  call comm_summation(ff1x,ff2x,lg_num(1)*ng_num(2)*ng_num(3),nproc_group_bound(1))

!$OMP parallel do private(kz,ky,kx)
  do kz = ng_sta(3),ng_end(3)
  do ky = ng_sta(2),ng_end(2)
  do kx = ng_sta(1),ng_end(1)
    ff1x(kx,ky,kz)=sum(eGxc(kx,:)*ff2x(:,ky,kz))/dble(lg_num(1)*lg_num(2)*lg_num(3))
  end do
  end do
  end do

  call comm_summation(ff1x,ff2x,lg_num(1)*ng_num(2)*ng_num(3),nproc_group_bound(1))

  bLx=2.d0*Pi/(Hgs(1)*dble(lg_num(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg_num(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg_num(3)))

!$OMP parallel do private(iz,iy,ix)
  do iz=lg_sta(3),lg_end(3)
  do iy=ng_sta(2),ng_end(2)
  do ix=ng_sta(1),ng_end(1)
    ff1z(ix,iy,iz)=0.d0
  end do
  end do
  end do

!$OMP parallel do private(kz,ky,kx,n,Gx,Gy,Gz,G2,kkx,kky,kkz)
  do kz = ng_sta(3),ng_end(3)
  do ky = ng_sta(2),ng_end(2)
  do kx = ng_sta(1),ng_end(1)
    n=(kz-lg_sta(3))*lg_num(2)*lg_num(1)+(ky-lg_sta(2))*lg_num(1)+kx-lg_sta(1)+1
!    Gx=2.d0*Pi*kx/lg_num(1)
!    Gy=2.d0*Pi*ky/lg_num(2)
!    Gz=2.d0*Pi*kz/lg_num(3)
    kkx=kx-1-lg_num(1)*(1+sign(1,(kx-1-lg_num(1)/2)))/2
    kky=ky-1-lg_num(2)*(1+sign(1,(ky-1-lg_num(2)/2)))/2
    kkz=kz-1-lg_num(3)*(1+sign(1,(kz-1-lg_num(3)/2)))/2
    Gx=kkx*bLx
    Gy=kky*bLy
    Gz=kkz*bLz
    G2=Gx**2+Gy**2+Gz**2
    if(kx-1==0.and.ky-1==0.and.kz-1==0)then
      rhoe_G_tmp(n)=0.d0
      ff1z(kx,ky,kz)=0.d0
    else
      rhoe_G_tmp(n)=ff2x(kx,ky,kz)
      ff1z(kx,ky,kz)=4.d0*Pi/G2*ff2x(kx,ky,kz)
    end if
  end do
  end do
  end do

  if(iSCFRT==1)then
    call comm_summation(rhoe_G_tmp,rhoe_G,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
  else if(iSCFRT==2)then
    if(itt==1.or.mod(itt,itcalc_ene)==0)then
      call comm_summation(rhoe_G_tmp,rhoe_G,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
    end if
  end if
  call comm_summation(ff1z,ff2z,ng_num(1)*ng_num(2)*lg_num(3),nproc_group_bound(3))

!$OMP parallel do private(iz,ky,kx)
  do iz = ng_sta(3),ng_end(3)
  do ky = ng_sta(2),ng_end(2)
  do kx = ng_sta(1),ng_end(1)
    ff1y(kx,ky,iz)=sum(eGz(:,iz)*ff2z(kx,ky,:))
  end do
  end do
  end do
  call comm_summation(ff1y,ff2y,ng_num(1)*lg_num(2)*ng_num(3),nproc_group_bound(2))

!$OMP parallel do private(iz,iy,kx)
  do iz = ng_sta(3),ng_end(3)
  do iy = ng_sta(2),ng_end(2)
  do kx = ng_sta(1),ng_end(1)
    ff1(kx,iy,iz)=sum(eGy(:,iy)*ff2y(kx,:,iz))
  end do
  end do
  end do
  call comm_summation(ff1,ff2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)

!$OMP parallel do private(iz,iy,ix)
  do iz = mg_sta(3),mg_end(3)
  do iy = mg_sta(2),mg_end(2)
  do ix = mg_sta(1),mg_end(1)
    tVh(ix,iy,iz)=sum(eGx(:,ix)*ff2(:,iy,iz))
  end do
  end do
  end do

  return
end subroutine Hartree_periodic
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120--------
