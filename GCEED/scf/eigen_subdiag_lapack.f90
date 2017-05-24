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

subroutine eigen_subdiag(Rmat,evec,iter,ier2)
use scf_data
implicit none

integer :: iter,ier2
real(8) :: Rmat(iter,iter)
real(8) :: evec(iter,iter)

character(1) :: JOBZ,UPLO
integer :: N
real(8) :: A(iter,iter)
integer :: LDA
real(8) :: W(iter)
real(8) :: WORK(3*iter-1)
integer :: LWORK

ier2=0

  JOBZ='V'
  UPLO='L'
  N=iter
  A=Rmat
  LDA=iter
  LWORK=3*iter-1
  call DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,ier2)
  evec=A

END subroutine eigen_subdiag
 
