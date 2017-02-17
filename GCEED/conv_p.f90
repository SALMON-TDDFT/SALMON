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

subroutine conv_p(iob,iob0)
use scf_data
implicit none
integer :: iob,iob0

if(ilsda==0)then
  iob0=iob
else if(ilsda==1)then
  if(iob<=MST(1))then
    iob0=iob
  else
    iob0=iob-MST(1)+MST0(1)
  end if
end if


end subroutine conv_p
