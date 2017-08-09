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
subroutine subdgemm(psi_mesh_box2,psibox,evec)
use scf_data
implicit none
integer :: N
real(8) :: psi_mesh_box2(ng_num(1)*ng_num(2)*ng_num(3),itotMST)
real(8) :: psibox(ng_num(1)*ng_num(2)*ng_num(3),itotMST)
real(8) :: evec(itotMST,itotMST)
real(8) :: ALPHA,BETA

N=ng_num(1)*ng_num(2)*ng_num(3)
ALPHA=1.d0
BETA=0.d0

call DGEMM('N','N',N,itotMST,itotMST,ALPHA, &
            psibox,N,evec,itotMST,BETA,psi_mesh_box2,N)

end subroutine subdgemm
