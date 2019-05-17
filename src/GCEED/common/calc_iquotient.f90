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
subroutine calc_iquotient(iob,nproc_ob,itotMST,iquotient)
implicit none
integer :: iob,nproc_ob,itotMST
integer :: iquotient

if(mod(iob*nproc_ob,itotMST)==0)then
  iquotient=iob*nproc_ob/itotMST-1
else
  iquotient=iob*nproc_ob/itotMST
end if

end subroutine calc_iquotient
