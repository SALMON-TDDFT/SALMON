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
subroutine copyVlocal(matbox12,ibox,is)
use scf_data
implicit none
integer :: ibox,is
real(8) :: matbox12(ista_Mxin_s(1,ibox):iend_Mxin_s(1,ibox),     &
                    ista_Mxin_s(2,ibox):iend_Mxin_s(2,ibox),     &
                    ista_Mxin_s(3,ibox):iend_Mxin_s(3,ibox))

  Vlocal( ista_Mxin_s(1,ibox):iend_Mxin_s(1,ibox),     &
          ista_Mxin_s(2,ibox):iend_Mxin_s(2,ibox),     &
          ista_Mxin_s(3,ibox):iend_Mxin_s(3,ibox),is) = &
           matbox12(ista_Mxin_s(1,ibox):iend_Mxin_s(1,ibox),     &
                    ista_Mxin_s(2,ibox):iend_Mxin_s(2,ibox),     &
                    ista_Mxin_s(3,ibox):iend_Mxin_s(3,ibox))

return

end subroutine copyVlocal
