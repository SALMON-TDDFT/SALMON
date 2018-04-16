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
subroutine initA(Ntime)
use scf_data
implicit none
integer :: Ntime
integer :: t_max

if(IC_rt==0)then
  t_max=Ntime
else
  t_max=Ntime+Miter_rt
end if

allocate( curr(3,0:t_max) )
allocate( sumcurr(3,0:t_max) )
allocate( A_ext(3,0:t_max+1) )
allocate( A_ind(3,0:t_max+1) )
allocate( A_tot(3,0:t_max+1) )
allocate( E_ext(3,0:t_max) )
allocate( E_ind(3,0:t_max) )
allocate( E_tot(3,0:t_max) )
allocate( rE_ind(3,0:t_max) )
curr=0.d0
sumcurr=0.d0
A_ind=0.d0
E_ind=0.d0
rE_ind=0.d0
A_ext(1:2,:)=0.d0
A_ext(3,:)=Fst
A_ind(:,0:1)=0.d0
A_tot(:,0:1)=A_ext(:,0:1)+A_ind(:,0:1)

E_ext=0.d0
E_ind(:,0)=0.d0
E_tot(:,0)=0.d0

end subroutine initA
