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

subroutine structure_opt
use scf_data
implicit none
integer :: iatom
real(8) :: alpha

alpha=1.d-1

do iatom=1,MI
  if(istopt_a(iatom)==1)then
    Rion(1:3,iatom)=Rion(1:3,iatom)+alpha*rforce(1:3,iatom)
  end if
end do

end subroutine structure_opt
