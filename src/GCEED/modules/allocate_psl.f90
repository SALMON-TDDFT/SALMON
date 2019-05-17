!
!  Copyright 2017-2019 SALMON developers
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
module allocate_psl_sub

use scf_data
use read_pslfile_sub

integer :: Nr

real(8), allocatable :: rhopp(:,:)
real(8), allocatable :: vpp(:,:,:)
real(8), allocatable :: uppr(:,:,:)
real(8), allocatable :: uVnl(:,:,:)

real(8), allocatable :: rad_psl(:,:)

real(8), allocatable :: ur(:,:)

real(8),allocatable :: rho_core(:,:,:)

real(8),allocatable :: vloctbl(:,:)
real(8),allocatable :: Gx(:),Gy(:),Gz(:)
real(8),allocatable :: dVloc_G(:,:)
real(8),allocatable :: dVloc_G_tmp(:,:)
complex(8),allocatable :: rhoion_G(:),Vion_G(:)
complex(8),allocatable :: rhoion_G_tmp(:)
complex(8),allocatable :: Vion_G_tmp(:)

complex(8),allocatable :: uVpsibox_c(:,:),uVpsibox2_c(:,:)
complex(8),allocatable :: uVpsibox1_j(:,:,:,:,:),uVpsibox2_j(:,:,:,:,:)

integer :: nGzero

contains
!==================================================================================================

subroutine allocate_psl
use scf_data
implicit none
real(8) :: r
integer :: Nr0

select case(iperiodic)
case(0)
  r=sqrt((dble(lg_num(1))*Hgs(1))**2  &
        +(dble(lg_num(2))*Hgs(2))**2  &
        +(dble(lg_num(3))*Hgs(3))**2)
case(3)
  r=sqrt((dble(lg_num(1))*Hgs(1))**2  &
        +(dble(lg_num(2))*Hgs(2))**2  &
        +(dble(lg_num(3))*Hgs(3))**2)+maxval(Rps(:))
end select
Nr0=r/rmin_step+1
if(Nr0<maxval(Mr(:)))then
  Nr=maxval(Mr(:))+1
else
  Nr=Nr0
end if

allocate(Mps_all(1:MI))
allocate(Mps(1:MI))

maxMps=int(4.d0/3.d0*Pi*(rmaxRps+4.d0*maxval(Hgs(:)))**3/Hvol)
Mlmps=(maxval(Mlps)+1)**2

allocate(Jxyz_all(1:3,1:maxMps,1:MI))
allocate(Jxyz_tmp1(1:3,1:maxMps,1:MI))
allocate(Jxyz_tmp2(1:3,1:maxMps,1:MI))
if(iperiodic==3)then
  allocate(Jxxyyzz_tmp1(1:3,1:maxMps,1:MI))
  allocate(Jxxyyzz_tmp2(1:3,1:maxMps,1:MI))
end if
allocate(uV_all(maxMps,Mlmps,MI), uVu(Mlmps,MI))

allocate(rhopp(0:Nr,MKI))
allocate(vpp(0:Nr,0:Nlps,MKI))
allocate(uppr(0:Nr,0:Nlps,MKI))
allocate(uVnl(0:Nr,0:Nlps,MKI))
allocate(rad_psl(0:Nr,MKI))

allocate(ur(maxMps,Nlmps))

allocate(rho_core(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

allocate(Jxyz(3,maxMps,MI))
allocate(uV(maxMps,Mlmps,MI))
allocate(numatom_ps(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

if(iperiodic==3)then
  allocate(vloctbl(Nr,MKI))
  allocate(Gx(lg_num(1)*lg_num(2)*lg_num(3)))
  allocate(Gy(lg_num(1)*lg_num(2)*lg_num(3)))
  allocate(Gz(lg_num(1)*lg_num(2)*lg_num(3)))
  allocate(dVloc_G(lg_num(1)*lg_num(2)*lg_num(3),MKI))
  allocate(dVloc_G_tmp(lg_num(1)*lg_num(2)*lg_num(3),MKI))
  allocate(rhoion_G(lg_num(1)*lg_num(2)*lg_num(3)))
  allocate(rhoion_G_tmp(lg_num(1)*lg_num(2)*lg_num(3)))
  allocate(Vion_G(lg_num(1)*lg_num(2)*lg_num(3)))
  allocate(Vion_G_tmp(lg_num(1)*lg_num(2)*lg_num(3)))
end if

allocate(Jxxyyzz_all(1:3,1:maxMps,1:MI))
allocate(Jxxyyzz(1:3,1:maxMps,1:MI))

if(iflag_ps==1) then
  allocate (uVpsibox_c(1:maxlm,1:MI))
  allocate (uVpsibox2_c(1:maxlm,1:MI))
  uVpsibox_c=0.d0
  uVpsibox2_c=0.d0
  if(iperiodic==3.and.iSCFRT==2)then
    allocate (uVpsibox1_j(4,1:maxlm,1:MI,1:iobnum,k_sta:k_end))
    allocate (uVpsibox2_j(4,1:maxlm,1:MI,1:iobnum,k_sta:k_end))
    uVpsibox1_j=0.d0
    uVpsibox2_j=0.d0
  end if
end if

end subroutine allocate_psl
!==================================================================================================
subroutine deallocate_psl

deallocate(Jxyz_all,Jxyz_tmp1,Jxyz_tmp2)
if(iperiodic==3) deallocate(Jxxyyzz_tmp1,Jxxyyzz_tmp2)
deallocate(Mps_all)
deallocate(Mps)
deallocate(uV_all,uVu)

deallocate(rhopp,vpp,uppr,uVnl)
deallocate(ur)
deallocate(rad_psl)

deallocate(rho_core)

deallocate(Jxyz,uV)
deallocate(numatom_ps)

if(iperiodic==3)then
  deallocate(vloctbl)
  deallocate(Gx,Gy,Gz)
  deallocate(dVloc_G,rhoion_G,Vion_G)
  deallocate(dVloc_G_tmp,rhoion_G_tmp,Vion_G_tmp)
end if
deallocate(Jxxyyzz_all)
deallocate(Jxxyyzz)

if(iflag_ps==1) then
  deallocate (uVpsibox_c,uVpsibox2_c)
  if(iperiodic==3.and.iSCFRT==2) deallocate (uVpsibox1_j,uVpsibox2_j)
end if

end subroutine deallocate_psl
!==================================================================================================
end module allocate_psl_sub
