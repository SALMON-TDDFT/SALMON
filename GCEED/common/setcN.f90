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
subroutine setcN
use scf_data
implicit none

cNmat(0,1)=-2.d0
cNmat(1,1)=1.d0

cNmat(0,2)=-5.d0/2.d0
cNmat(1,2)=4.d0/3.d0
cNmat(2,2)=-1.d0/12.d0

cNmat(0,3)=-49.d0/18.d0
cNmat(1,3)=3.d0/2.d0
cNmat(2,3)=-3.d0/20.d0
cNmat(3,3)=1.d0/90.d0

cNmat(0,5)=-5269.d0/1800.d0
cNmat(1,5)=5.d0/3.d0
cNmat(2,5)=-5.d0/21.d0
cNmat(3,5)=5.d0/126.d0
cNmat(4,5)=-5.d0/1008.d0
cNmat(5,5)=1.d0/3150.d0

cNmat(0,4)=-205.d0/72.d0
cNmat(1,4)=8.d0/5.d0
cNmat(2,4)=-1.d0/5.d0
cNmat(3,4)=8.d0/315.d0
cNmat(4,4)=-1.d0/560.d0

cNmat(0,6)=-5369.d0/1800.d0
cNmat(1,6)=12.d0/7.d0
cNmat(2,6)=-15.d0/56.d0
cNmat(3,6)=10.d0/189.d0
cNmat(4,6)=-1.d0/112.d0
cNmat(5,6)=2.d0/1925.d0
cNmat(6,6)=-1.d0/16632.d0

cNmat(0,7)=-266681.d0/88200.d0
cNmat(1,7)=7.d0/4.d0
cNmat(2,7)=-7.d0/24.d0
cNmat(3,7)=7.d0/108.d0
cNmat(4,7)=-7.d0/528.d0
cNmat(5,7)=7.d0/3300.d0
cNmat(6,7)=-7.d0/30888.d0
cNmat(7,7)=1.d0/84084.d0

cNmat(0,8)=-1077749.d0/352800.d0
cNmat(1,8)=16.d0/9.d0
cNmat(2,8)=-14.d0/45.d0
cNmat(3,8)=112.d0/1485.d0
cNmat(4,8)=-7.d0/396.d0
cNmat(5,8)=112.d0/32175.d0
cNmat(6,8)=-2.d0/3861.d0
cNmat(7,8)=16.d0/315315.d0
cNmat(8,8)=-1.d0/411840.d0

cNmat(0,9)=-9778141.d0/3175200.d0
cNmat(1,9)=9.d0/5.d0
cNmat(2,9)=-18.d0/55.d0
cNmat(3,9)=14.d0/165.d0
cNmat(4,9)=-63.d0/2860.d0
cNmat(5,9)=18.d0/3575.d0
cNmat(6,9)=-2.d0/2145.d0
cNmat(7,9)=9.d0/70070.d0
cNmat(8,9)=-9.d0/777920.d0
cNmat(9,9)=1.d0/1969110.d0

cNmat(0,10)=-1968329.d0/635040.d0
cNmat(1,10)=20.d0/11.d0
cNmat(2,10)=-15.d0/44.d0
cNmat(3,10)=40.d0/429.d0
cNmat(4,10)=-15.d0/572.d0
cNmat(5,10)=24.d0/3575.d0
cNmat(6,10)=-5.d0/3432.d0
cNmat(7,10)=30.d0/119119.d0
cNmat(8,10)=-5.d0/155584.d0
cNmat(9,10)=10.d0/3741309.d0
cNmat(10,10)=-1.d0/9237800.d0

cNmat(0,11)=-239437889.d0/76839840.d0
cNmat(1,11)=11.d0/6.d0
cNmat(2,11)=-55.d0/156.d0
cNmat(3,11)=55.d0/546.d0
cNmat(4,11)=-11.d0/364.d0
cNmat(5,11)=11.d0/1300.d0
cNmat(6,11)=-11.d0/5304.d0
cNmat(7,11)=55.d0/129948.d0
cNmat(8,11)=-55.d0/806208.d0
cNmat(9,11)=11.d0/1360476.d0
cNmat(10,11)=-11.d0/17635800.d0
cNmat(11,11)=1.d0/42678636.d0

cNmat(0,12)=-240505109.d0/76839840.d0
cNmat(1,12)=24.d0/13.d0
cNmat(2,12)=-33.d0/91.d0
cNmat(3,12)=88.d0/819.d0
cNmat(4,12)=-99.d0/2912.d0
cNmat(5,12)=396.d0/38675.d0
cNmat(6,12)=-11.d0/3978.d0
cNmat(7,12)=132.d0/205751.d0
cNmat(8,12)=-33.d0/268736.d0
cNmat(9,12)=44.d0/2380833.d0
cNmat(10,12)=-3.d0/1469650.d0
cNmat(11,12)=12.d0/81800719.d0
cNmat(12,12)=-1.d0/194699232.d0


end subroutine setcN
