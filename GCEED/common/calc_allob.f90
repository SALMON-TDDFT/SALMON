! Copyright 2017 Katsuyuki Nobusada, Masashi Noda, Kazuya Ishimura, Kenji Iida, Maiku Yamaguchi
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

subroutine calc_allob(iob,iob_allob)
use scf_data
use new_world_sub
implicit none
integer :: iob,iob_allob

if(ilsda==0.or.nproc_ob==1)then
  if(iparaway_ob==1)then
    iob_allob=newrank_comm_grid*itotMST/nproc_ob+iob
  else if(iparaway_ob==2)then
    iob_allob=(iob-1)*nproc_ob+newrank_comm_grid+1
  end if
else
  if(iparaway_ob==1)then
    if(newrank_comm_spin<nproc_ob_spin(1))then
      iob_allob=newrank_comm_grid*MST(1)/nproc_ob_spin(1)+iob
    else
      iob_allob=newrank_comm_grid*MST(2)/nproc_ob_spin(2)+iob+MST(1)
    end if
  else if(iparaway_ob==2)then
    if(newrank_comm_spin<nproc_ob_spin(1))then
      iob_allob=(iob-1)*nproc_ob_spin(1)+newrank_comm_grid+1
    else
      iob_allob=(iob-1)*nproc_ob_spin(2)+newrank_comm_grid+1+MST(1)
    end if
  end if
end if

end subroutine calc_allob
