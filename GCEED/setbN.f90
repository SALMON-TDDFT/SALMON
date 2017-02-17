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

subroutine setbN
use scf_data

bNmat(1,1)=1.d0/2.d0

bNmat(1,2)=2.d0/3.d0
bNmat(2,2)=-1.d0/12.d0

bNmat(1,3)=3.d0/4.d0
bNmat(2,3)=-3.d0/20.d0
bNmat(3,3)=1.d0/60.d0

bNmat(1,4)=4.d0/5.d0
bNmat(2,4)=-1.d0/5.d0
bNmat(3,4)=4.d0/105.d0
bNmat(4,4)=-1.d0/280.d0

end subroutine setbN

