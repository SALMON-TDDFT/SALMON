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
subroutine calc_myob(iob,iob_myob)
use salmon_parallel, only: nproc_id_spin
use scf_data
use new_world_sub
implicit none
integer :: iob,iob_myob
integer :: iquotient,iob_min

if(ilsda==0.or.nproc_ob==1)then
  if(iparaway_ob==1)then
    call calc_iquotient(iob,nproc_ob,itotMST,iquotient)
    iob_min=itotMST*iquotient/nproc_ob
    iob_myob=iob-iob_min
  else if(iparaway_ob==2)then
    iob_myob=(iob-1)/nproc_ob+1
  end if
else
  if(iparaway_ob==1)then
    if(nproc_id_spin<nproc_ob_spin(1))then
      call calc_iquotient(iob,nproc_ob_spin(1),MST(1),iquotient)
      iob_min=MST(1)*iquotient/nproc_ob_spin(1)
      iob_myob=iob-iob_min
    else
      call calc_iquotient(iob-MST(1),nproc_ob_spin(2),MST(2),iquotient)
      iob_min=MST(2)*iquotient/nproc_ob_spin(2)
      iob_myob=iob-MST(1)-iob_min
    end if
  else if(iparaway_ob==2)then
    if(nproc_id_spin<nproc_ob_spin(1))then
      iob_myob=(iob-1)/nproc_ob_spin(1)+1
    else
      iob_myob=(iob-MST(1)-1)/nproc_ob_spin(2)+1
    end if
  end if
end if

end subroutine calc_myob
