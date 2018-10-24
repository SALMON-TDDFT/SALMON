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
subroutine eigen_subdiag_periodic(Rmat,evec,iter,ier2)
  use scf_data
  implicit none
  character :: JOBZ, UPLO
  integer :: LWORK
  integer :: iter,ier2
  real(8),allocatable :: RWORK(:)
  real(8) :: W(iter)
  complex(8) :: Rmat(iter,iter)
  complex(8),allocatable :: WORK(:)
  complex(8) :: evec(iter,iter)
  
  ier2=0
  
  JOBZ='V'
  UPLO='U'
  
  LWORK=2*iter-1
  allocate(WORK(LWORK))
  allocate(RWORK(3*iter-2))
  
  call ZHEEV(JOBZ,UPLO,iter,Rmat,iter,W,WORK,LWORK,RWORK,ier2)
  
  evec(:,:)=Rmat(:,:)
  
  deallocate(WORK,RWORK)
  
end subroutine eigen_subdiag_periodic
