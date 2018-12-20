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
MODULE allocate_mat_sub

use inputoutput, only: iperiodic
use scf_data
implicit none

real(8), allocatable :: vecR(:,:,:,:)
real(8), allocatable :: vecR_tmp(:,:,:,:)

real(8), allocatable :: wk_s_h(:,:,:),wk2_s_h(:,:,:),lap_wk_s_h(:,:,:)
real(8), allocatable :: wkbound_h(:),wk2bound_h(:)

real(8), allocatable :: matbox_m(:,:,:),matbox_m2(:,:,:)
complex(8), allocatable :: cmatbox_m(:,:,:),cmatbox_m2(:,:,:)
real(8), allocatable :: matbox_l(:,:,:),matbox_l2(:,:,:)
complex(8), allocatable :: cmatbox_l(:,:,:),cmatbox_l2(:,:,:)

complex(8), allocatable :: zalpha2(:,:,:,:),zalpha3(:,:,:,:)

real(8),allocatable :: rgrad_wk(:,:,:,:,:,:)

complex(8),allocatable :: cgrad_wk(:,:,:,:,:,:)

real(8), allocatable :: rho_tmp(:,:,:)
real(8), allocatable :: rho_s_tmp(:,:,:,:)
real(8), allocatable :: vxc_tmp(:,:,:)
real(8), allocatable :: vxc_s_tmp(:,:,:,:)
real(8), allocatable :: eexc_tmp(:,:,:)
real(8), allocatable :: exc_dummy(:,:,:)
real(8), allocatable :: exc_dummy2(:,:,:,:)
real(8), allocatable :: exc_dummy3(:,:,:,:)

complex(8),allocatable :: rhoe_G(:)
complex(8),allocatable :: rhoe_G_tmp(:)

! hartree routine for periodic system
real(8),allocatable :: trho2z(:,:,:)
real(8),allocatable :: trho3z(:,:,:)
complex(8),allocatable :: eGx(:,:)
complex(8),allocatable :: eGxc(:,:)
complex(8),allocatable :: eGy(:,:)
complex(8),allocatable :: eGyc(:,:)
complex(8),allocatable :: eGz(:,:)
complex(8),allocatable :: eGzc(:,:)
complex(8),allocatable :: ff1(:,:,:)
complex(8),allocatable :: ff2(:,:,:)
complex(8),allocatable :: ff1x(:,:,:)
complex(8),allocatable :: ff2x(:,:,:)
complex(8),allocatable :: ff1y(:,:,:)
complex(8),allocatable :: ff2y(:,:,:)
complex(8),allocatable :: ff1z(:,:,:)
complex(8),allocatable :: ff2z(:,:,:)

! FFTE routine
complex(8),allocatable :: A_FFTE(:,:,:), B_FFTE(:,:,:)
real(8),allocatable :: coef_poisson(:,:,:)
real(8),allocatable :: coef_poisson_FFTE(:,:,:)
real(8),allocatable :: A_FFTE_copy(:,:,:), A_FFTE_copy2(:,:,:)

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE allocate_mat

implicit none

allocate (vecR(3,lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )

allocate (vecR_tmp(3,lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )

allocate (matbox_m(mg_sta(1):mg_end(1),    &
             mg_sta(2):mg_end(2),      &
             mg_sta(3):mg_end(3)) )

allocate (matbox_m2(mg_sta(1):mg_end(1),    &
             mg_sta(2):mg_end(2),      &
             mg_sta(3):mg_end(3)) )
allocate (cmatbox_m(mg_sta(1):mg_end(1),    &
             mg_sta(2):mg_end(2),      &
             mg_sta(3):mg_end(3)) )
allocate (cmatbox_m2(mg_sta(1):mg_end(1),    &
             mg_sta(2):mg_end(2),      &
             mg_sta(3):mg_end(3)) )

allocate (matbox_l(lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )
allocate (matbox_l2(lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )
allocate (cmatbox_l(lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )
allocate (cmatbox_l2(lg_sta(1):lg_end(1),    &
             lg_sta(2):lg_end(2),      &
             lg_sta(3):lg_end(3)) )

allocate (wk_s_h(ng_sta(1)-Ndh:ng_end(1)+Ndh,   &
             ng_sta(2)-Ndh:ng_end(2)+Ndh,   &
             ng_sta(3)-Ndh:ng_end(3)+Ndh))
allocate (wk2_s_h(ng_sta(1):ng_end(1),   &
              ng_sta(2):ng_end(2),   &
              ng_sta(3):ng_end(3)))
allocate (lap_wk_s_h(ng_sta(1):ng_end(1),   &
                 ng_sta(2):ng_end(2),   &
                 ng_sta(3):ng_end(3)))

allocate (wkbound_h(lg_num(1)*lg_num(2)*lg_num(3)/minval(lg_num(1:3))*6*Ndh) )
allocate (wk2bound_h(lg_num(1)*lg_num(2)*lg_num(3)/minval(lg_num(1:3))*6*Ndh) )

if(icalcforce==1)then
  allocate(rforce(3,MI))
end if

if(iSCFRT==1.and.icalcforce==1)then
  select case(iperiodic)
  case(0)
    allocate(rgrad_wk(mg_sta(1):mg_end(1)+1,   &
                      mg_sta(2):mg_end(2),     &
                      mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end,3))
  case(3)
    allocate(cgrad_wk(mg_sta(1):mg_end(1)+1,   &
                      mg_sta(2):mg_end(2),     &
                      mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end,3))
  end select
else if(iSCFRT==2.and.icalcforce==1)then
  allocate(cgrad_wk(mg_sta(1):mg_end(1)+1,   &
                    mg_sta(2):mg_end(2),     &
                    mg_sta(3):mg_end(3),1:iobnum,k_sta:k_end,3))
end if

if(iSCFRT==2)then
  iwk_size=2
  call make_iwksta_iwkend

  if(iflag_fourier_omega==1)then
    allocate(zalpha2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3),num_fourier_omega))
    allocate(zalpha3(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3),num_fourier_omega))
  end if

  if(iflag_md==1)then
    allocate(dRion(3,MI,-1:1))
    allocate(Rion_eq(3,MI))
    dRion(:,:,:)=0.d0
    Rion_eq(:,:)=Rion(:,:)
  end if

end if

if(iSCFRT==1.and.iperiodic==0)then
  allocate (rxk_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
  allocate (rhxk_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
  allocate (rgk_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
  allocate (rpk_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
end if

if(iSCFRT==1.and.iperiodic==3)then
  allocate (zxk_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
  allocate (zhxk_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
  allocate (zgk_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
  allocate (zpk_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
  allocate (zpko_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
  allocate (zhtpsi_ob(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),1:iobnum))
end if

allocate (rho_tmp(ng_num(1), ng_num(2), ng_num(3)))
allocate (rho_s_tmp(ng_num(1), ng_num(2), ng_num(3), 2))
allocate (vxc_tmp(ng_num(1), ng_num(2), ng_num(3)))
allocate (vxc_s_tmp(ng_num(1), ng_num(2), ng_num(3), 2))
allocate (eexc_tmp(ng_num(1), ng_num(2), ng_num(3)))
allocate (exc_dummy(ng_num(1), ng_num(2), ng_num(3)))
allocate (exc_dummy2(ng_num(1), ng_num(2), ng_num(3),2))
allocate (exc_dummy3(ng_num(1), ng_num(2), ng_num(3),3))

select case(iperiodic)
case(3)
  if(iSCFRT==2.and.iflag_hartree==4)then
    allocate(rhoe_G(lg_num(1)*lg_num(2)/NPUY*lg_num(3)/NPUZ))
  else
    allocate(rhoe_G(lg_num(1)*lg_num(2)*lg_num(3)))
    allocate(rhoe_G_tmp(lg_num(1)*lg_num(2)*lg_num(3)))
  end if

  if(iflag_hartree==2)then
    allocate(trho2z(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),lg_sta(3):lg_end(3)))
    allocate(trho3z(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),lg_sta(3):lg_end(3)))
    allocate(eGx(lg_sta(1):lg_end(1),lg_sta(1):lg_end(1)))
    allocate(eGxc(lg_sta(1):lg_end(1),lg_sta(1):lg_end(1)))
    allocate(eGy(lg_sta(2):lg_end(2),lg_sta(2):lg_end(2)))
    allocate(eGyc(lg_sta(2):lg_end(2),lg_sta(2):lg_end(2)))
    allocate(eGz(lg_sta(3):lg_end(3),lg_sta(3):lg_end(3)))
    allocate(eGzc(lg_sta(3):lg_end(3),lg_sta(3):lg_end(3)))
    allocate(ff1(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
    allocate(ff2(lg_sta(1):lg_end(1),lg_sta(2):lg_end(2),lg_sta(3):lg_end(3)))
    allocate(ff1x(lg_sta(1):lg_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
    allocate(ff2x(lg_sta(1):lg_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
    allocate(ff1y(ng_sta(1):ng_end(1),lg_sta(2):lg_end(2),ng_sta(3):ng_end(3)))
    allocate(ff2y(ng_sta(1):ng_end(1),lg_sta(2):lg_end(2),ng_sta(3):ng_end(3)))
    allocate(ff1z(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),lg_sta(3):lg_end(3)))
    allocate(ff2z(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),lg_sta(3):lg_end(3)))
    allocate(coef_poisson(ng_sta(1):ng_end(1),ng_sta(2):ng_end(2),ng_sta(3):ng_end(3)))
  else if(iflag_hartree==4)then
    allocate(A_FFTE(lg_num(1),lg_num(2)/NPUY,lg_num(3)/NPUZ))
    allocate(B_FFTE(lg_num(1),lg_num(2)/NPUY,lg_num(3)/NPUZ))
    allocate(coef_poisson(lg_num(1),lg_num(2)/NPUY,lg_num(3)/NPUZ))
    allocate(A_FFTE_copy(lg_num(1),lg_num(2),lg_num(3)/NPUZ))
    allocate(A_FFTE_copy2(lg_num(1),lg_num(2),lg_num(3)/NPUZ))
  end if
end select

allocate(icoo1d(3,lg_num(1)*lg_num(2)*lg_num(3)))

END SUBROUTINE allocate_mat

!======================================================================

END MODULE allocate_mat_sub
