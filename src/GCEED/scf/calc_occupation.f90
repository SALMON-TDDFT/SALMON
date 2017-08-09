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
subroutine calc_occupation
use scf_data
implicit none
integer :: is_sta,is_end
integer :: pstart(2),pend(2)

if(ilsda==0)then
  is_sta=1
  is_end=1
else
  is_sta=1
  is_end=2
  pstart(1)=1
  pend(1)=MST(1)
  pstart(2)=MST(1)+1
  pend(2)=itotMST
end if

rocc(1:itotMST,1)=0.d0
if(ilsda==0)then
  rocc(1:ifMST(1),1)=2.d0
else
  rocc(1:ifMST(1),1)=1.d0
  rocc(MST(1)+1:MST(1)+ifMST(2),1)=1.d0
end if

end subroutine calc_occupation

