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
MODULE deallocate_mat_sub

use scf_data
use allocate_mat_sub
implicit none

CONTAINS

!=======================================================================
!=======================================================================

SUBROUTINE deallocate_mat

implicit none

deallocate (vecR)
deallocate (vecR_tmp)
deallocate (matbox_m,matbox_m2)
deallocate (matbox_l,matbox_l2)
deallocate (cmatbox_m,cmatbox_m2)
deallocate (cmatbox_l,cmatbox_l2)

deallocate (wk_s_h,wk2_s_h,lap_wk_s_h)
deallocate (wkbound_h,wk2bound_h)

if(icalcforce==1) then
  deallocate(rforce)
end if

if(iSCFRT==2)then
  if(iflag_fourier_omega==1)then
    deallocate(zalpha2,zalpha3)
  end if
  if(iflag_MD==1)then
    deallocate(dRion,Rion_eq)
  end if
end if

if(iSCFRT==1.and.icalcforce==1)then
  deallocate(rgrad_wk)
else if(iSCFRT==2.and.icalcforce==1)then
  deallocate(cgrad_wk)
end if

if(iSCFRT==1.and.iperiodic==0)then
  deallocate(rxk_ob,rhxk_ob,rgk_ob,rpk_ob)
end if

if(iSCFRT==1.and.iperiodic==3)then
  deallocate(zxk_ob,zhxk_ob,zgk_ob,zpk_ob)
  deallocate(zpko_ob,zhtpsi_ob)
end if

deallocate (rho_tmp)
deallocate (rho_s_tmp)
deallocate (vxc_tmp)
deallocate (vxc_s_tmp)
deallocate (eexc_tmp)
deallocate (exc_dummy)
deallocate (exc_dummy2)
deallocate (exc_dummy3)

deallocate (icoo1d)

END SUBROUTINE deallocate_mat

!======================================================================

END MODULE deallocate_mat_sub
