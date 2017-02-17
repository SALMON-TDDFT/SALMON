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

!=====================================================================
subroutine calc_invA_p(Amat,ix)
use scf_data
implicit none

integer :: ii,jj,ix
integer :: M,N,LDA
integer :: INFO
integer,allocatable :: IPIV(:)
complex(8) :: Amat(max_lg_num,max_lg_num)
complex(8),allocatable :: WORK(:)
complex(8),allocatable :: Amat_lib(:,:)
integer :: LWORK

allocate(Amat_lib(lg_num(ix),lg_num(ix)))

do jj=1,lg_num(ix)
do ii=1,lg_num(ix)
  Amat_lib(ii,jj)=Amat(ii,jj)
end do
end do

M=lg_num(ix)
N=lg_num(ix)
LDA=M
allocate(IPIV(N))

call ZGETRF( M, N, Amat_lib, LDA, IPIV, INFO )

allocate(WORK(N))
LWORK=N

call ZGETRI( N, Amat_lib, LDA, IPIV, WORK, LWORK, INFO )

do jj=1,lg_num(ix)
do ii=1,lg_num(ix)
  Amat(ii,jj)=Amat_lib(ii,jj)
end do
end do

deallocate(IPIV)
deallocate(Amat_lib)

end subroutine calc_invA_p

!=====================================================================
subroutine calc_invA_psl(Amat,iatom)
use scf_data
implicit none

integer :: iatom,ii,jj
integer :: M,N,LDA
integer :: INFO
integer,allocatable :: IPIV(:)
complex(8) :: Amat(maxMps,maxMps)
complex(8),allocatable :: WORK(:)
complex(8),allocatable :: Amat_lib(:,:)
integer :: LWORK

allocate(Amat_lib(Mps(iatom),Mps(iatom)))

do jj=1,Mps(iatom)
do ii=1,Mps(iatom)
  Amat_lib(ii,jj)=Amat(ii,jj)
end do
end do

M=Mps(iatom)
N=Mps(iatom)
LDA=M
allocate(IPIV(N))

call ZGETRF( M, N, Amat_lib, LDA, IPIV, INFO )

allocate(WORK(N))
LWORK=N

call ZGETRI( N, Amat_lib, LDA, IPIV, WORK, LWORK, INFO )

do jj=1,Mps(iatom)
do ii=1,Mps(iatom)
  Amat(ii,jj)=Amat_lib(ii,jj)
end do
end do

deallocate(IPIV)
deallocate(Amat_lib)

end subroutine calc_invA_psl

