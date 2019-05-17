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
subroutine calc_iroot(iob,iroot)
use salmon_parallel, only: nproc_id_spin
use scf_data
use new_world_sub
implicit none
integer :: iob,iroot

if(ilsda==0.or.nproc_ob==1)then
  if(iparaway_ob==1)then
    call calc_iquotient(iob,nproc_ob,itotMST,iroot)
  else if(iparaway_ob==2)then
    iroot=mod(iob-1,nproc_ob)
  end if
else
  if(iparaway_ob==1)then
    if(nproc_id_spin<nproc_ob_spin(1))then
      call calc_iquotient(iob,nproc_ob_spin(1),MST(1),iroot)
    else
      call calc_iquotient(iob-MST(1),nproc_ob_spin(2),MST(2),iroot)
    end if
  else if(iparaway_ob==2)then
    if(nproc_id_spin<nproc_ob_spin(1))then
      iroot=mod(iob-1,nproc_ob_spin(1))
    else
      iroot=mod(iob-1-MST(1),nproc_ob_spin(2))
    end if
  end if 
end if

end subroutine calc_iroot
