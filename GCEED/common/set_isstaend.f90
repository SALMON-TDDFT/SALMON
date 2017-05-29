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
subroutine set_isstaend(is_sta,is_end)
use scf_data
use new_world_sub
implicit none
integer :: is_sta,is_end

if(ilsda==0)then
  is_sta=1
  is_end=1
else
  if(nproc_ob==1)then
    is_sta=1
    is_end=2
  else
    if(newrank_comm_spin<nproc_ob_spin(1))then
      is_sta=1
      is_end=1
    else
      is_sta=2
      is_end=2
    end if
  end if
end if

return

end subroutine set_isstaend
