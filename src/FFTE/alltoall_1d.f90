!
!  Copyright 2018 SALMON developers
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
!======================================================================
subroutine alltoall_1d_y(vec_a,vec_b,nn,npuy)
  use salmon_communication, only: comm_alltoall
  use salmon_parallel, only: nproc_group_icommy
  implicit none
  complex(8),intent(in)   :: vec_a(nn)
  complex(8),intent(out)  :: vec_b(nn)
  integer,intent(in)      :: nn
  integer,intent(in)      :: npuy

  call comm_alltoall(vec_a,vec_b,nproc_group_icommy,nn/npuy)

end subroutine alltoall_1d_y
!======================================================================
subroutine alltoall_1d_z(vec_a,vec_b,nn,npuz)
  use salmon_communication, only: comm_alltoall
  use salmon_parallel, only: nproc_group_icommz
  implicit none
  complex(8),intent(in)   :: vec_a(nn)
  complex(8),intent(out)  :: vec_b(nn)
  integer,intent(in)      :: nn
  integer,intent(in)      :: npuz

  call comm_alltoall(vec_a,vec_b,nproc_group_icommz,nn/npuz)

end subroutine alltoall_1d_z
