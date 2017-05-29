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
MODULE eigen_sub

INTERFACE eigenval
  MODULE PROCEDURE R_eigenval, C_eigenval
END INTERFACE

INTERFACE gen_eigen
  MODULE PROCEDURE R_gen_eigen, C_gen_eigen
END INTERFACE

CONTAINS

!=======================================================================
!======================================================= RMM-DIIS_lapack
!==================================================== eigenvalue problem
! This is a routine to solve eigenvalue problem.
subroutine R_eigenval(Smat, eval, iter)

implicit none
integer :: ii,iter
real(8) :: Smat(iter,iter),eval(iter)
character :: JOBVL*1, JOBVR*1
integer :: INFO, LDA, LDB, LDVL, LDVR, LWORK, N
real(8),allocatable :: A( :, : ), ALPHAI( : ), ALPHAR( : )
real(8),allocatable :: B( :, : ), BETA( : ), VL( :, : )
real(8),allocatable :: VR( :, : ), WORK( : )

JOBVL='N'
JOBVR='V'
N=iter
LDA=iter
LDB=iter
     
LDVL=iter
LDVR=iter

LWORK=8*iter

allocate ( A( LDA, N ), ALPHAI( N ), ALPHAR( N ) )
allocate ( B( LDB, N ), BETA( N ), VL( LDVL, N ) )
allocate ( VR( LDVR, N ), WORK( LWORK ) )

A=Smat

B=0.d0
do ii=1,iter
  B(ii,ii)=1.d0
end do

call DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,      &
       BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

! check whether eigenvalues are real.

do ii=1,iter
  if(ALPHAI(ii) > 0.05d0)then
    write(*,*) "========== DIIS error =========="
    write(*,*) "One of eigenvalue is imaginary."
    write(*,*) "================================"
    stop
  end if
end do

eval=ALPHAR

deallocate (A,ALPHAI,ALPHAR,B,BETA,VL,VR,WORK)

return

end subroutine R_eigenval

! This is a routine to solve generalised eigenvalue problem.

SUBROUTINE R_gen_eigen(Rmat,Smat,alpha,betav,evec,iter,ier2)

implicit none

integer :: iter,ii,ier2
real(8) :: Rmat(iter,iter)
real(8) :: Smat(iter,iter)
real(8) :: alpha(iter),betav(iter)
real(8) :: evec(iter,iter)
character :: JOBVL*1, JOBVR*1
integer :: INFO, LDA, LDB, LDVL, LDVR, LWORK, N
real(8),allocatable :: A( :, : ), ALPHAI( : ), ALPHAR( : )
real(8),allocatable :: B( :, : ), BETA( : ), VL( :, : )
real(8),allocatable :: VR( :, : ), WORK( : )

ier2=0

JOBVL='N'
JOBVR='V'
N=iter
LDA=iter
LDB=iter
     
LDVL=iter
LDVR=iter

LWORK=8*iter

allocate ( A( LDA, N ), ALPHAI( N ), ALPHAR( N ) )
allocate ( B( LDB, N ), BETA( N ), VL( LDVL, N ) )
allocate ( VR( LDVR, N ), WORK( LWORK ) )

A=Rmat
B=Smat

call DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,      &
       BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

! check whether eigenvalues and eigenvectors take real value.
do ii=1,iter
  if(ALPHAI(ii) > 1.d-8)then
    write(*,*) "========== DIIS error =========="
    write(*,*) "One of eigenvector is imaginary."
    write(*,*) "================================"
    stop
  end if
end do

alpha=ALPHAR
betav=BETA
evec=VR

deallocate (A,ALPHAI,ALPHAR,B,BETA,VL,VR,WORK)

return

end subroutine R_gen_eigen

!=======================================================================
!======================================================= RMM-DIIS_lapack
!==================================================== eigenvalue problem
! This is a routine to solve eigenvalue problem.
subroutine C_eigenval(Smat, eval, iter)

implicit none
integer :: ii,iter
complex(8) :: Smat(iter,iter)
real(8) :: eval(iter)
character :: JOBVL*1, JOBVR*1
integer :: INFO, LDA, LDB, LDVL, LDVR, LWORK, N
complex(8),allocatable :: A( :, : ), ALPHA( : )
complex(8),allocatable :: B( :, : ), BETA( : ), VL( :, : )
complex(8),allocatable :: VR( :, : ), WORK( : )

JOBVL='N'
JOBVR='V'
N=iter
LDA=iter
LDB=iter
     
LDVL=iter
LDVR=iter

LWORK=8*iter

allocate ( A( LDA, N ), ALPHA( N ) )
allocate ( B( LDB, N ), BETA( N ), VL( LDVL, N ) )
allocate ( VR( LDVR, N ), WORK( LWORK ) )

A=Smat

B=0.d0
do ii=1,iter
  B(ii,ii)=1.d0
end do

call ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA,      &
       BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

eval=real(ALPHA)

deallocate (A,ALPHA,B,BETA,VL,VR,WORK)

return

end subroutine C_eigenval

! This is a routine to solve generalised eigenvalue problem.

SUBROUTINE C_gen_eigen(Rmat,Smat,alpha,betav,evec,iter,ier2)

implicit none

integer :: iter,ier2
complex(8) :: Rmat(iter,iter)
complex(8) :: Smat(iter,iter)
complex(8) :: alpha(iter),betav(iter)
complex(8) :: evec(iter,iter)
character :: JOBVL*1, JOBVR*1
integer :: INFO, LDA, LDB, LDVL, LDVR, LWORK, N
complex(8),allocatable :: A( :, : ), ALPHA2( : )
complex(8),allocatable :: B( :, : ), BETA( : ), VL( :, : )
complex(8),allocatable :: VR( :, : ), WORK( : )

ier2=0

JOBVL='N'
JOBVR='V'
N=iter
LDA=iter
LDB=iter
     
LDVL=iter
LDVR=iter

LWORK=8*iter

allocate ( A( LDA, N ), ALPHA2( N ) )
allocate ( B( LDB, N ), BETA( N ), VL( LDVL, N ) )
allocate ( VR( LDVR, N ), WORK( LWORK ) )

A=Rmat
B=Smat

call ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA2,      &
       BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )

alpha=ALPHA2
betav=BETA
evec=VR

deallocate (A,ALPHA2,B,BETA,VL,VR,WORK)

return

end subroutine C_gen_eigen

END MODULE eigen_sub
