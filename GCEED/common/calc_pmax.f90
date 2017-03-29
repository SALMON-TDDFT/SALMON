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

subroutine calc_pmax(iobmax)
use scf_data
use new_world_sub
implicit none
integer :: iobmax,iss

if(iSCFRT==1)then
  iobmax=iobnum
else if(iSCFRT==2)then
  if(ilsda==0)then
    call calc_iobnum(ifMST(1),nproc_ob,newrank_comm_grid,iobmax,nproc_ob,iparaway_ob)
  else if(ilsda==1)then
    iobmax=iobnum
  end if
end if

end subroutine calc_pmax
