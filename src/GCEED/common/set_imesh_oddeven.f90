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
subroutine set_imesh_oddeven(itmg)
  use scf_data
  !$use omp_lib
  implicit none
  integer :: itmg
  integer :: jj

  do jj=1,3
    if(mod(int(rLsize(jj,itmg)/Hgs(jj)+1.d-12),2)==1)then
      imesh_oddeven(jj)=1
    else
      imesh_oddeven(jj)=2
    end if
  end do 
  
end subroutine set_imesh_oddeven
