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
subroutine taylor_coe
use scf_data
implicit none
integer :: mm,nn
complex(8),parameter :: zi=(0.d0,1.d0)

do nn=1,N_hamil
   zc(nn)=(-zi*dt)**nn
   do mm=1,nn
      zc(nn)=zc(nn)/mm
   end do
end do

end subroutine taylor_coe
