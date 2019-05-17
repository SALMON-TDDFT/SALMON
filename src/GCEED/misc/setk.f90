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
subroutine setk(k_sta, k_end, k_num, num_kpoints_rd, nproc_k, newrank_comm_orbitalgrid)

implicit none
integer :: k_sta,k_end,k_num
integer :: nproc_k
integer :: num_kpoints_rd
integer :: newrank_comm_orbitalgrid

k_sta=(newrank_comm_orbitalgrid)*num_kpoints_rd/nproc_k+1
k_end=(newrank_comm_orbitalgrid+1)*num_kpoints_rd/nproc_k
k_num=(newrank_comm_orbitalgrid+1)*num_kpoints_rd/nproc_k-(newrank_comm_orbitalgrid)*num_kpoints_rd/nproc_k

end subroutine setk
