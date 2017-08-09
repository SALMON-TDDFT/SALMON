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
subroutine calc_allob(iob,iob_allob)
use salmon_parallel, only: nproc_id_grid, nproc_id_spin
use scf_data
use new_world_sub
implicit none
integer :: iob,iob_allob

if(ilsda==0.or.nproc_ob==1)then
  if(iparaway_ob==1)then
    iob_allob=nproc_id_grid*itotMST/nproc_ob+iob
  else if(iparaway_ob==2)then
    iob_allob=(iob-1)*nproc_ob+nproc_id_grid+1
  end if
else
  if(iparaway_ob==1)then
    if(nproc_id_spin<nproc_ob_spin(1))then
      iob_allob=nproc_id_grid*MST(1)/nproc_ob_spin(1)+iob
    else
      iob_allob=nproc_id_grid*MST(2)/nproc_ob_spin(2)+iob+MST(1)
    end if
  else if(iparaway_ob==2)then
    if(nproc_id_spin<nproc_ob_spin(1))then
      iob_allob=(iob-1)*nproc_ob_spin(1)+nproc_id_grid+1
    else
      iob_allob=(iob-1)*nproc_ob_spin(2)+nproc_id_grid+1+MST(1)
    end if
  end if
end if

end subroutine calc_allob
