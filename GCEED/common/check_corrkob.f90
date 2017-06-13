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
subroutine check_corrkob(iob,icorr_p)
use salmon_parallel, only: nproc_id_grid, nproc_id_spin
use scf_data
use new_world_sub
implicit none
integer :: iob,icorr_p
integer :: iquotient

if(ilsda==0.or.nproc_ob==1)then
  if(iparaway_ob==1)then
    call calc_iquotient(iob,nproc_ob,itotMST,iquotient)
    if(nproc_id_grid==iquotient)then
      icorr_p=1
    else
      icorr_p=0
    end if
  else if(iparaway_ob==2)then
    if(nproc_id_grid==mod(iob-1,nproc_ob))then
      icorr_p=1
    else
      icorr_p=0
    end if
  end if
else
  if(iparaway_ob==1)then
    if(nproc_id_spin<nproc_ob_spin(1))then
      call calc_iquotient(iob,nproc_ob_spin(1),MST(1),iquotient)
      if(iob<=MST(1).and.nproc_id_grid==iquotient)then
        icorr_p=1
      else
        icorr_p=0
      end if
    else
      call calc_iquotient(iob-MST(1),nproc_ob_spin(2),MST(2),iquotient)
      if(iob>=MST(1)+1.and.nproc_id_grid==iquotient)then
        icorr_p=1
      else
        icorr_p=0
      end if
    end if
  else if(iparaway_ob==2)then
    if(nproc_id_spin<nproc_ob_spin(1))then
      if(iob<=MST(1).and.nproc_id_grid==mod(iob-1,nproc_ob_spin(1)))then
        icorr_p=1
      else
        icorr_p=0
      end if
    else
      if(iob>=MST(1)+1.and.nproc_id_grid==mod(iob-1-MST(1),nproc_ob_spin(2)))then
        icorr_p=1
      else
        icorr_p=0
      end if
    end if
  end if
end if

end subroutine check_corrkob
