!
!  Copyright 2018 SALMON developers
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
subroutine calcVpsl_periodic
  use salmon_parallel, only: nproc_group_global, nproc_size_global, nproc_id_global, nproc_group_grid
  use salmon_parallel, only: nproc_id_icommy, nproc_id_icommz
  use salmon_communication, only: comm_bcast, comm_summation, comm_is_root
  use scf_data
  use new_world_sub
  use allocate_psl_sub
  use allocate_mat_sub
  implicit none
  
  integer :: ii,ix,iy,iz,ak
  integer :: iix,iiy,iiz
  integer :: iatom
  real(8) :: x,y,z
  
  integer :: n
  real(8) :: bLx,bLy,bLz
  real(8) :: aLxyz
  integer :: NG_s,NG_e
  integer :: NG_l_s_para,NG_l_e_para
  integer :: numtmp
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8) :: G2sq,G2
  real(8) :: Gd
  real(8) :: Gr
  real(8) :: s
  real(8) :: r
  integer :: imax
  real(8) :: Vpsl_tmp(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))

!calculate reciprocal lattice vector
  bLx=2.d0*Pi/(Hgs(1)*dble(lg_num(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg_num(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg_num(3)))
  
 
  NG_s=1
  NG_e=lg_num(1)*lg_num(2)*lg_num(3)
  
  numtmp=(NG_e-NG_s+1)/nproc_size_global
  
  NG_l_s_para = nproc_id_global*numtmp+1
  NG_l_e_para = (nproc_id_global+1)*numtmp
  if(nproc_id_global==nproc_size_global-1) NG_l_e_para=NG_e
  
  nGzero=-1
  
  do ak=1,MKI
    do ii=1,Mr(ak)
      vloctbl(ii,ak)=vpp(ii,Lref(ak),ak)
    enddo
  end do
  
  n=0
  do iz=lg_sta(3),lg_end(3)
  do iy=lg_sta(2),lg_end(2)
  do ix=lg_sta(1),lg_end(1)
    n=n+1
    if((ix-1)**2+(iy-1)**2+(iz-1)**2 == 0) nGzero=n
    iix=ix-1-lg_num(1)*(1+sign(1,(ix-1-(lg_num(1)+1)/2)))/2
    iiy=iy-1-lg_num(2)*(1+sign(1,(iy-1-(lg_num(2)+1)/2)))/2
    iiz=iz-1-lg_num(3)*(1+sign(1,(iz-1-(lg_num(3)+1)/2)))/2
    Gx(n)=dble(iix)*bLx
    Gy(n)=dble(iiy)*bLy
    Gz(n)=dble(iiz)*bLz
  enddo
  enddo
  enddo

! local potential
  dVloc_G_tmp(:,:)=0.d0
  do ak=1,MKI
    imax=min(Mr(ak),Nr-1)
    do n=NG_l_s_para, NG_l_e_para
      G2sq=sqrt(Gx(n)**2+Gy(n)**2+Gz(n)**2)
      s=0.d0
      if (n == nGzero) then
        do ii=2,imax
          r=rad_psl(ii,ak)
          s=s+4*Pi*r**2*(vloctbl(ii,ak)+Zps(ak)/r)*(rad_psl(ii+1,ak)-rad_psl(ii,ak))
        enddo
      else
        do ii=2,imax
          r=rad_psl(ii,ak)
          s=s+4*Pi*r**2*sin(G2sq*r)/(G2sq*r)*(vloctbl(ii,ak)+Zps(ak)/r)*(rad_psl(ii+1,ak)-rad_psl(ii,ak))
        enddo
      endif
      dVloc_G_tmp(n,ak)=s
    enddo
  enddo
  call comm_summation(dVloc_G_tmp,dVloc_G,(NG_e-NG_s+1)*MKI,nproc_group_global)
 
  aLxyz=Hvol*dble(lg_num(1)*lg_num(2)*lg_num(3))
  rhoion_G_tmp=0.d0
  do iatom=1,MI
    do n=NG_l_s_para,NG_l_e_para
      rhoion_G_tmp(n)=rhoion_G_tmp(n)+Zps(Kion(iatom))/aLxyz*exp(-zI*(Gx(n)*Rion(1,iatom)+Gy(n)*Rion(2,iatom)+Gz(n)*Rion(3,iatom)))
    enddo
  enddo
  call comm_summation(rhoion_G_tmp,rhoion_G,NG_e-NG_s+1,nproc_group_global)

  Vion_G_tmp=0.d0
  do iatom=1,MI
    ak=Kion(iatom)
    do n=NG_l_s_para,NG_l_e_para
      G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
      Gd=Gx(n)*Rion(1,iatom)+Gy(n)*Rion(2,iatom)+Gz(n)*Rion(3,iatom)
      Vion_G_tmp(n)=Vion_G_tmp(n)+dVloc_G(n,ak)*exp(-zI*Gd)/aLxyz
      if(n == nGzero) cycle
      Vion_G_tmp(n)=Vion_G_tmp(n)-4*Pi/G2*Zps(ak)*exp(-zI*Gd)/aLxyz
    enddo
  enddo
  call comm_summation(Vion_G_tmp,Vion_G,NG_e-NG_s+1,nproc_group_global)

  Vpsl_tmp=0.d0
  do n=NG_s,NG_e
  !$OMP parallel do private(iz,iy,ix,x,y,z,Gr)
    do iz=ng_sta(3),ng_end(3)
    do iy=ng_sta(2),ng_end(2)
    do ix=ng_sta(1),ng_end(1)
      x=gridcoo(ix,1)
      y=gridcoo(iy,2)
      z=gridcoo(iz,3)
      Gr=Gx(n)*x+Gy(n)*y+Gz(n)*z
      Vpsl_tmp(ix,iy,iz)=Vpsl_tmp(ix,iy,iz)+Vion_G(n)*exp(zI*Gr)
    enddo
    enddo
    enddo
  enddo
  call comm_summation(Vpsl_tmp,Vpsl,mg_num(1)*mg_num(2)*mg_num(3),nproc_group_grid)

  return

end subroutine calcVpsl_periodic
