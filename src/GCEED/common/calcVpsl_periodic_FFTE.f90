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
subroutine calcVpsl_periodic_FFTE
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
  
  integer :: n
  real(8) :: bLx,bLy,bLz
  real(8) :: aLxyz
  integer :: NG_s,NG_e
  integer :: NG_l_s_para,NG_l_e_para
  integer :: numtmp
  complex(8),parameter :: zI=(0.d0,1.d0)
  real(8) :: G2sq,G2
  real(8) :: Gd
  real(8) :: s
  real(8) :: r
  integer :: imax
  integer :: iy_sta,iy_end,iz_sta,iz_end
  integer :: i,iix2,iiy2,iiz2
  

!calculate reciprocal lattice vector
  bLx=2.d0*Pi/(Hgs(1)*dble(lg_num(1)))
  bLy=2.d0*Pi/(Hgs(2)*dble(lg_num(2)))
  bLz=2.d0*Pi/(Hgs(3)*dble(lg_num(3)))

  iz_sta=1
  iz_end=lg_num(3)/NPUZ
  iy_sta=1
  iy_end=lg_num(2)/NPUY
 
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

  do iz=1,lg_num(3)/NPUZ
  do iy=1,lg_num(2)/NPUY
  do ix=1,lg_num(1)
    n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
    iix2=ix-1+lg_sta(1)
    iiy2=iy-1+nproc_id_icommy*lg_num(2)/NPUY+lg_sta(2)
    iiz2=iz-1+nproc_id_icommz*lg_num(3)/NPUZ+lg_sta(3)
    if(ix==1.and.iy==1.and.iz==1.and.nproc_id_icommz==0.and.nproc_id_icommy==0) nGzero=n
    iix=ix-1-lg_num(1)*(1+sign(1,(iix2-1-(lg_num(1)+1)/2)))/2
    iiy=iy-1+nproc_id_icommy*lg_num(2)/NPUY-lg_num(2)*(1+sign(1,(iiy2-1-(lg_num(2)+1)/2)))/2
    iiz=iz-1+nproc_id_icommz*lg_num(3)/NPUZ-lg_num(3)*(1+sign(1,(iiz2-1-(lg_num(3)+1)/2)))/2
    Gx(n)=dble(iix)*bLx
    Gy(n)=dble(iiy)*bLy
    Gz(n)=dble(iiz)*bLz
  enddo
  enddo
  enddo

  dVloc_G(:,:)=0.d0
  do ak=1,MKI
    imax=min(Mr(ak),Nr-1)
    do iz=1,lg_num(3)/NPUZ
    do iy=1,lg_num(2)/NPUY
    do ix=1,lg_num(1)
      n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
      G2sq=sqrt(Gx(n)**2+Gy(n)**2+Gz(n)**2)
      s=0.d0
      if (n == nGzero) then
        do i=2,imax
          r=rad_psl(i,ak)
          s=s+4*Pi*r**2*(vloctbl(i,ak)+Zps(ak)/r)*(rad_psl(i+1,ak)-rad_psl(i,ak))
        enddo
      else
        do i=2,imax
          r=rad_psl(i,ak)
          s=s+4*Pi*r**2*sin(G2sq*r)/(G2sq*r)*(vloctbl(i,ak)+Zps(ak)/r)*(rad_psl(i+1,ak)-rad_psl(i,ak))
        enddo
      endif
      dVloc_G(n,ak)=s
    enddo
    enddo
    enddo
  enddo
 
  aLxyz=Hvol*dble(lg_num(1)*lg_num(2)*lg_num(3))
  rhoion_G=0.d0
  do iatom=1,MI
    do iz=1,lg_num(3)/NPUZ
    do iy=1,lg_num(2)/NPUY
    do ix=1,lg_num(1)
      n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
      rhoion_G(n)=rhoion_G(n)+Zps(Kion(iatom))/aLxyz*exp(-zI*(Gx(n)*Rion(1,iatom)+Gy(n)*Rion(2,iatom)+Gz(n)*Rion(3,iatom)))
    enddo
    enddo
    enddo
  enddo

  Vion_G=0.d0
  do iatom=1,MI
    ak=Kion(iatom)
    do iz=1,lg_num(3)/NPUZ
    do iy=1,lg_num(2)/NPUY
    do ix=1,lg_num(1)
      n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
      G2=Gx(n)**2+Gy(n)**2+Gz(n)**2
      Gd=Gx(n)*Rion(1,iatom)+Gy(n)*Rion(2,iatom)+Gz(n)*Rion(3,iatom)
      Vion_G(n)=Vion_G(n)+dVloc_G(n,ak)*exp(-zI*Gd)/aLxyz
      if(n == nGzero) cycle
      Vion_G(n)=Vion_G(n)-4*Pi/G2*Zps(ak)*exp(-zI*Gd)/aLxyz
    enddo
    enddo
    enddo
  enddo

  CALL PZFFT3DV_MOD(A_FFTE,B_FFTE,lg_num(1),lg_num(2),lg_num(3),NPUY,NPUZ,0) 

  do iz=1,lg_num(3)/NPUZ
  do iy=1,lg_num(2)/NPUY
  do ix=1,lg_num(1)
    n=(iz-1)*lg_num(2)/NPUY*lg_num(1)+(iy-1)*lg_num(1)+ix
    B_FFTE(ix,iy,iz)=Vion_G(n)*dble(lg_num(1)*lg_num(2)*lg_num(3))
  enddo
  enddo
  enddo

  CALL PZFFT3DV_MOD(B_FFTE,A_FFTE,lg_num(1),lg_num(2),lg_num(3),NPUY,NPUZ,1)

  if(icheck_ascorder==1)then
!$OMP parallel do
    do iz = lg_sta(3),lg_end(3)
    do iy = lg_sta(2),lg_end(2)
    do ix = lg_sta(1),lg_end(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do
    if(NPUW==1)then
!$OMP parallel do private(iiz,iiy)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg_num(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg_num(2)/NPUY
!            Vpsl(1:lg_end(1),iiy,iiz)=A_FFTE(1:lg_end(1),iy,iz)
          matbox_l(1:lg_end(1),iiy,iiz)=A_FFTE(1:lg_end(1),iy,iz)
        end do
      end do
    else
!$OMP parallel do private(iiz,iiy,ix)
      do iz=iz_sta,iz_end
        iiz=iz+nproc_id_icommz*lg_num(3)/NPUZ
        do iy=iy_sta,iy_end
          iiy=iy+nproc_id_icommy*lg_num(2)/NPUY
          do iix=ng_sta(1),ng_end(1)
            ix=iix-lg_sta(1)+1
            matbox_l(iix,iiy,iiz)=A_FFTE(ix,iy,iz)
          end do
        end do
      end do
    end if
    call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
!$OMP parallel do
    do iz = mg_sta(3),mg_end(3)
    do iy = mg_sta(2),mg_end(2)
    do ix = mg_sta(1),mg_end(1)
      Vpsl(ix,iy,iz)=matbox_l2(ix,iy,iz) 
    end do
    end do
    end do
  else
!$OMP parallel do
    do iz = lg_sta(3),lg_end(3)
    do iy = lg_sta(2),lg_end(2)
    do ix = lg_sta(1),lg_end(1)
      matbox_l(ix,iy,iz)=0.d0
    end do
    end do
    end do

!$OMP parallel do private(iiz,iiy)
    do iz=iz_sta,iz_end
      iiz=iz+nproc_id_icommz*lg_num(3)/NPUZ
      do iy=iy_sta,iy_end
        iiy=iy+nproc_id_icommy*lg_num(2)/NPUY
        matbox_l(1:lg_end(1),iiy,iiz)=A_FFTE(1:lg_end(1),iy,iz)
      end do
    end do
  
    call comm_summation(matbox_l,matbox_l2,lg_num(1)*lg_num(2)*lg_num(3),nproc_group_global)
   
!$OMP parallel do
    do iz = mg_sta(3),mg_end(3)
    do iy = mg_sta(2),mg_end(2)
    do ix = mg_sta(1),mg_end(1)
      Vpsl(ix,iy,iz)=matbox_l2(ix,iy,iz)
    end do
    end do
    end do
  end if

  return

end subroutine calcVpsl_periodic_FFTE
