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

r=sqrt((dble(lg_num(1))*Hgs(1))**2  &
      +(dble(lg_num(2))*Hgs(2))**2  &
      +(dble(lg_num(3))*Hgs(3))**2)
Nr0=r/rmin_step+1
if(Nr0<maxval(Mr(:)))then
  Nr=maxval(Mr(:))+1
else
  Nr=Nr0
end if

allocate(Mps(1:MI))

maxMps=int(4.d0/3.d0*Pi*(rmaxRps+4.d0*maxval(Hgs(:)))**3/Hvol)
Mlmps=(maxval(Mlps)+1)**2

allocate(Jxyz(1:3,1:maxMps,1:MI))
allocate(Jxyz_tmp1(1:3,1:maxMps,1:MI))
allocate(Jxyz_tmp2(1:3,1:maxMps,1:MI))
allocate(uV(maxMps,Mlmps,MI), uVu(Mlmps,MI))

allocate(rhopp(0:Nr,MKI))
allocate(vpp(0:Nr,0:Nlps,MKI))
allocate(uppr(0:Nr,0:Nlps,MKI))
allocate(uVnl(0:Nr,0:Nlps,MKI))
allocate(rad_psl(0:Nr,MKI))

allocate(ur(maxMps,Nlmps))

allocate(rho_core(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

allocate(jMps_l(maxMps,MI))
allocate(max_jMps_l(MI))

allocate(jMps_l_s(maxMps,MI))
allocate(max_jMps_l_s(MI))

allocate(Jxyz2nd(3,maxMps,MI))
allocate(uV2nd(maxMps,Mlmps,MI))
allocate(numatom_ps(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3)))

allocate(Jxxyyzz(1:3,1:maxMps,1:MI))
allocate(Jxxyyzz2nd(1:3,1:maxMps,1:MI))

if(iflag_ps==1) then
  allocate (uVpsibox_c(1:maxlm,1:MI))
  allocate (uVpsibox2_c(1:maxlm,1:MI))
  uVpsibox_c=0.d0
  uVpsibox2_c=0.d0
end if

end subroutine allocate_psl
!==================================================================================================
subroutine deallocate_psl

deallocate(Jxyz,Jxyz_tmp1,Jxyz_tmp2)
deallocate(Mps)
deallocate(uV,uVu)

deallocate(rhopp,vpp,uppr,uVnl)
deallocate(ur)
deallocate(rad_psl)

deallocate(rho_core)

deallocate(jMps_l,max_jMps_l)
deallocate(jMps_l_s,max_jMps_l_s)

deallocate(Jxyz2nd,uV2nd)
deallocate(numatom_ps)

deallocate(Jxxyyzz)
deallocate(Jxxyyzz2nd)

if(iflag_ps==1) then
  deallocate (uVpsibox_c,uVpsibox2_c)
end if

end subroutine deallocate_psl
!==================================================================================================
end module allocate_psl_sub
