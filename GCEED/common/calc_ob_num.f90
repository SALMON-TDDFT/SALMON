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
subroutine calc_iobnum(tMST,tproc,trank,tiobnum,nproc_ob,iparaway_ob)
implicit none
integer :: tMST,tproc,trank,tiobnum
integer :: nproc_ob,iparaway_ob

if(iparaway_ob==1)then
  tiobnum=(trank+1)*tMST/nproc_ob-trank*tMST/nproc_ob
else if(iparaway_ob==2)then
  if(mod(tMST,tproc)==0)then
    tiobnum=tMST/tproc
  else
    if(trank<mod(tMST,tproc))then
      tiobnum=tMST/tproc+1
    else
      tiobnum=tMST/tproc
    end if
  end if
end if

end subroutine calc_iobnum
